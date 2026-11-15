//! A local HTTP server that plays the GitHub Releases API and serves release files, for tests
//! (`FIELDCRAFT_UPDATE_HOST=http://127.0.0.1:PORT`). It answers each request from a table of
//! paths and counts the requests, so a test can also check that nothing was asked.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, PoisonError};

type Routes = Arc<Mutex<HashMap<String, (u16, Vec<u8>)>>>;

/// The server; it runs until the process ends.
pub struct Server {
    port: u16,
    routes: Routes,
    requests: Arc<AtomicUsize>,
}

impl Server {
    /// Start on a free local port. `None` when no port can be bound.
    pub fn start() -> Option<Server> {
        let listener = TcpListener::bind("127.0.0.1:0").ok()?;
        let port = listener.local_addr().ok()?.port();
        let routes: Routes = Arc::default();
        let requests = Arc::new(AtomicUsize::new(0));
        let (r, n) = (routes.clone(), requests.clone());
        std::thread::Builder::new()
            .name("fieldcraft-update mock".into())
            .spawn(move || {
                for stream in listener.incoming().flatten() {
                    n.fetch_add(1, Ordering::SeqCst);
                    let _ = answer(stream, &r);
                }
            })
            .ok()?;
        Some(Server { port, routes, requests })
    }

    /// Its base URL, for `FIELDCRAFT_UPDATE_HOST`.
    pub fn url(&self) -> String {
        format!("http://127.0.0.1:{}", self.port)
    }

    /// Answer `path` with `status` and `body`.
    pub fn route(&self, path: &str, status: u16, body: Vec<u8>) {
        self.routes.lock().unwrap_or_else(PoisonError::into_inner).insert(path.to_string(), (status, body));
    }

    /// Answer the API's `releases/latest` (`which` = `latest`) or `releases/tags/<which>` with
    /// `release`.
    pub fn release(&self, which: &str, release: &serde_json::Value) {
        let path = if which == "latest" {
            format!("/repos/{}/releases/latest", crate::REPO)
        } else {
            format!("/repos/{}/releases/tags/{which}", crate::REPO)
        };
        self.route(&path, 200, release.to_string().into_bytes());
    }

    /// Answer the API's release list (`releases`) with `list`.
    pub fn releases(&self, list: &serde_json::Value) {
        self.route(&format!("/repos/{}/releases", crate::REPO), 200, list.to_string().into_bytes());
    }

    /// How many requests it has had.
    pub fn requests(&self) -> usize {
        self.requests.load(Ordering::SeqCst)
    }
}

fn answer(stream: TcpStream, routes: &Routes) -> std::io::Result<()> {
    let mut reader = BufReader::new(stream.try_clone()?);
    let mut line = String::new();
    reader.read_line(&mut line)?;
    // The headers, up to the blank line.
    let mut header = String::new();
    while reader.read_line(&mut header)? > 2 {
        header.clear();
    }
    let path = line.split_whitespace().nth(1).unwrap_or("/").to_string();
    let (status, body) =
        routes.lock().unwrap_or_else(PoisonError::into_inner).get(&path).cloned().unwrap_or((404, b"{\"message\":\"Not Found\"}".to_vec()));
    let mut stream = stream;
    write!(stream, "HTTP/1.1 {status} X\r\nContent-Length: {}\r\nContent-Type: application/octet-stream\r\nConnection: close\r\n\r\n", body.len())?;
    stream.write_all(&body)?;
    stream.flush()
}
