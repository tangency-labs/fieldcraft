//! Bundled, user-selected skills. Names and resources are a compile-time allowlist; fetching
//! one never reaches a backend, opens a file or runs a modelling command.

use serde_json::{Value, json};

struct Prompt {
    name: &'static str,
    title: &'static str,
    description: &'static str,
    text: &'static str,
    resources: &'static [Resource],
}

struct Resource {
    uri: &'static str,
    name: &'static str,
    mime: &'static str,
    text: &'static str,
}

const REFERENCES: &[Resource] = &[
    Resource {
        uri: "fieldcraft://skills/fieldcraft-modelling/references/checks.md",
        name: "Modelling verification exercises",
        mime: "text/markdown",
        text: include_str!("../../../skills/fieldcraft-modelling/references/checks.md"),
    },
    Resource {
        uri: "fieldcraft://skills/fieldcraft-modelling/references/axial-rod.json",
        name: "Original axial rod recipe",
        mime: "application/json",
        text: include_str!("../../../skills/fieldcraft-modelling/references/axial-rod.json"),
    },
    Resource {
        uri: "fieldcraft://skills/fieldcraft-modelling/references/steady-conductor.json",
        name: "Original steady conductor recipe",
        mime: "application/json",
        text: include_str!("../../../skills/fieldcraft-modelling/references/steady-conductor.json"),
    },
];

const PROMPTS: &[Prompt] = &[
    Prompt {
        name: "fieldcraft-meshing",
        title: "FieldCraft meshing",
        description: "Plan, inspect and verify a mesh with the bundled meshing skill.",
        text: include_str!("../../../skills/fieldcraft-meshing/SKILL.md"),
        resources: &[],
    },
    Prompt {
        name: "fieldcraft-modelling",
        title: "FieldCraft modelling",
        description: "Choose and check a structural, heat or fluid model with the bundled modelling skill.",
        text: include_str!("../../../skills/fieldcraft-modelling/SKILL.md"),
        resources: REFERENCES,
    },
];

fn object(params: &Value) -> Result<&serde_json::Map<String, Value>, String> {
    params.as_object().ok_or_else(|| "prompt parameters must be an object".into())
}

/// This small immutable catalogue is returned in one page. An absent or null cursor requests
/// that page; a non-null continuation cursor cannot name a page and is refused.
pub(crate) fn list(params: &Value) -> Result<Value, String> {
    if !params.is_null() && object(params)?.get("cursor").is_some_and(|cursor| !cursor.is_null()) {
        return Err("the bundled prompt catalogue has no continuation cursor".into());
    }
    Ok(json!({"prompts": PROMPTS.iter().map(|p| json!({
        "name": p.name, "title": p.title, "description": p.description, "arguments": [],
    })).collect::<Vec<_>>()}))
}

pub(crate) fn get(params: &Value) -> Result<Value, String> {
    let params = object(params)?;
    let name = params.get("name").and_then(Value::as_str).ok_or_else(|| "missing prompt `name` (a string)".to_string())?;
    let prompt = PROMPTS
        .iter()
        .find(|p| p.name == name)
        .ok_or_else(|| format!("unknown prompt `{}`; use prompts/list for the bundled skills", name.chars().take(200).collect::<String>()))?;
    if let Some(arguments) = params.get("arguments")
        && !arguments.as_object().is_some_and(|a| a.is_empty())
    {
        return Err("bundled skill prompts take no arguments; omit `arguments` or send an empty object".into());
    }
    let mut messages = vec![json!({"role": "user", "content": {"type": "text", "text": prompt.text}})];
    // Keep the skill payload byte-identical. Linked, fixed references are separate MCP
    // embedded resources, also discoverable/readable through resources/list and resources/read.
    for resource in prompt.resources {
        messages.push(json!({"role": "user", "content": {"type": "resource", "resource": resource.content()}}));
    }
    Ok(json!({"description": prompt.description, "messages": messages}))
}

impl Resource {
    fn content(&self) -> Value {
        json!({"uri": self.uri, "mimeType": self.mime, "text": self.text})
    }
}

pub(crate) fn resources() -> Vec<Value> {
    REFERENCES.iter().map(|r| json!({"uri": r.uri, "name": r.name, "mimeType": r.mime})).collect()
}

pub(crate) fn read(uri: &str) -> Option<Value> {
    REFERENCES.iter().find(|r| r.uri == uri).map(|r| json!({"contents": [r.content()]}))
}
