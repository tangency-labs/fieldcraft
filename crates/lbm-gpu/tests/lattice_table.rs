//! The shader's D3Q19 direction table (`C` and `c_of` in lbm.wgsl) inside a loop over the
//! directions, the pattern NVIDIA's Vulkan driver (GB10, 580.159.03) miscompiled when the
//! directions came from a `switch` or from a table of vec3<i32>: the wall-function bounce's loop takes ½ τ_x d_x off each
//! population of a wall cell whose d = c_q − n is one step along x (n the wall's normal, ±y). The
//! table and `c_of` are cut from the solver's own shader, so this checks the code the solver runs.
//! Skipped without a GPU adapter unless `FIELDCRAFT_REQUIRE_GPU=1`; `FIELDCRAFT_GPU` picks it.

use std::borrow::Cow;

use fieldcraft_gpu::select::{self, Policy, gpu_required};
use fieldcraft_lbm::lattice::C;

const LOOP: &str = r#"
@group(0) @binding(0) var<storage, read> normals: array<i32>;
@group(0) @binding(1) var<storage, read_write> f: array<f32>;

@compute @workgroup_size(64)
fn shift(@builtin(global_invocation_id) id: vec3<u32>) {
    let k = id.x;
    if k >= arrayLength(&normals) { return; }
    let n = vec3<i32>(0, normals[k], 0);
    let tau = 1e-3 * f32(k + 1u);
    let cells = arrayLength(&normals);
    for (var q = 1u; q < 19u; q = q + 1u) {
        let d = c_of(q) - n;
        if dot(d, d) != 1 || d.x == 0 { continue; }
        f[q * cells + k] = f[q * cells + k] - 0.5 * tau * f32(d.x);
    }
}
"#;

/// `const C` through the end of `c_of` in the solver's shader.
fn table() -> &'static str {
    let shader = include_str!("../src/lbm.wgsl");
    let start = shader.find("const C:").expect("lbm.wgsl defines the direction table `C`");
    let end = start + shader[start..].find("fn opp(").expect("lbm.wgsl defines `opp` after `c_of`");
    &shader[start..end]
}

#[test]
fn the_direction_table_gives_every_direction_in_a_loop() {
    let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
    let adapter = match select::adapter(policy) {
        Ok(Some((adapter, _))) => adapter,
        Ok(None) | Err(_) if !gpu_required() => {
            println!("skipping: no GPU adapter");
            return;
        }
        other => panic!("a GPU is required: {:?}", other.err()),
    };
    println!("GPU: {}", adapter.get_info().name);
    let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default())).expect("a device");
    // Wall cells 0..8 with the normal +y (a bottom wall), 8..16 with −y (a top wall).
    let normals: Vec<i32> = (0..16).map(|k| if k < 8 { 1 } else { -1 }).collect();
    let cells = normals.len();
    let pops: Vec<f32> = (0..19 * cells).map(|i| 0.01 * (i % 97) as f32).collect();

    let source = format!("{}{LOOP}", table());
    let module = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: None, source: wgpu::ShaderSource::Wgsl(Cow::Owned(source)) });
    let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: None,
        layout: None,
        module: &module,
        entry_point: Some("shift"),
        compilation_options: wgpu::PipelineCompilationOptions::default(),
        cache: None,
    });
    let storage = |bytes: Vec<u8>| {
        let buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: bytes.len() as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&buffer, 0, &bytes);
        buffer
    };
    let normals_buffer = storage(normals.iter().flat_map(|v| v.to_le_bytes()).collect());
    let f = storage(pops.iter().flat_map(|v| v.to_le_bytes()).collect());
    let bind = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: None,
        layout: &pipeline.get_bind_group_layout(0),
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: normals_buffer.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: f.as_entire_binding() },
        ],
    });
    let read = device.create_buffer(&wgpu::BufferDescriptor {
        label: None,
        size: f.size(),
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor::default());
        pass.set_pipeline(&pipeline);
        pass.set_bind_group(0, &bind, &[]);
        pass.dispatch_workgroups(1, 1, 1);
    }
    encoder.copy_buffer_to_buffer(&f, 0, &read, 0, f.size());
    queue.submit([encoder.finish()]);
    read.slice(..).map_async(wgpu::MapMode::Read, |_| {});
    device.poll(wgpu::PollType::wait_indefinitely()).expect("the device finishes");
    let got: Vec<f32> = read.slice(..).get_mapped_range().expect("mapped").chunks(4).map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]])).collect();

    let mut want = pops.clone();
    for (k, &ny) in normals.iter().enumerate() {
        let tau = 1e-3 * (k + 1) as f32;
        for (q, c) in C.iter().enumerate().skip(1) {
            let d = [c[0], c[1] - ny, c[2]];
            if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] == 1 && d[0] != 0 {
                want[q * cells + k] -= 0.5 * tau * d[0] as f32;
            }
        }
    }
    let wrong: Vec<(usize, usize)> =
        (0..19).flat_map(|q| (0..cells).map(move |k| (q, k))).filter(|&(q, k)| got[q * cells + k] != want[q * cells + k]).collect();
    assert!(wrong.is_empty(), "populations (q, cell) wrong: {wrong:?}");
}
