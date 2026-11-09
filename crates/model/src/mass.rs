/// Mass matrix choice for a linear dynamic analysis.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MassKind {
    #[default]
    Consistent,
    Lumped,
}
