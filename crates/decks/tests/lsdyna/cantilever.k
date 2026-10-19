*KEYWORD long=y
$ Written by FieldCraft: implicit linear analysis, long format (20-column fields).
$ Units: mm-t-s-N-MPa (millimetre, tonne, second, newton, megapascal, degree Celsius, consistent).
*CONTROL_UNITS
$             length                time                mass                temp
                  mm                 sec            mtrc_ton                   c
*NODE
$                nid                   x                   y                   z
                   1                 0.0                 0.0                 0.0
                   2                 1.0                 0.0                 0.0
                   3                 2.0                 0.0                 0.0
                   4                 3.0                 0.0                 0.0
                   5                 4.0                 0.0                 0.0
                   6                 0.0                 1.0                 0.0
                   7                 1.0                 1.0                 0.0
                   8                 2.0                 1.0                 0.0
                   9                 3.0                 1.0                 0.0
                  10                 4.0                 1.0                 0.0
                  11                 0.0                 0.0                 1.0
                  12                 1.0                 0.0                 1.0
                  13                 2.0                 0.0                 1.0
                  14                 3.0                 0.0                 1.0
                  15                 4.0                 0.0                 1.0
                  16                 0.0                 1.0                 1.0
                  17                 1.0                 1.0                 1.0
                  18                 2.0                 1.0                 1.0
                  19                 3.0                 1.0                 1.0
                  20                 4.0                 1.0                 1.0
*MAT_ELASTIC_TITLE
Steel
$                mid                  ro                   e                  pr
                   1             7.85e-9            210000.0                 0.3
*SECTION_SOLID
$              secid              elform
                   1                   2
*PART
beam
$                pid               secid                 mid               eosid                hgid
                   1                   1                   1
*ELEMENT_SOLID
$                eid                 pid
$                 n1                  n2                  n3                  n4                  n5                  n6                  n7                  n8                  n9                 n10
                   1                   1
                   1                   2                   7                   6                  11                  12                  17                  16
                   2                   1
                   2                   3                   8                   7                  12                  13                  18                  17
                   3                   1
                   3                   4                   9                   8                  13                  14                  19                  18
                   4                   1
                   4                   5                  10                   9                  14                  15                  20                  19
*SET_NODE_LIST_TITLE
root
$                sid
                   2
                   1                   6                  11                  16
*SET_NODE_LIST_TITLE
tip
$                sid
                   3
                   5                  10                  15                  20
*TITLE
Bend
*CONTROL_IMPLICIT_GENERAL
$             imflag                 dt0
                   1                 1.0
*CONTROL_IMPLICIT_SOLUTION
$             nsolvr
                   1
*CONTROL_TERMINATION
$             endtim
                 1.0
*DEFINE_CURVE_TITLE
Bend ramp
$               lcid
                   1
$                 a1                  o1
                 0.0                 0.0
                 1.0                 1.0
*BOUNDARY_SPC_SET_ID
                   1clamp
$               nsid                 cid                dofx                dofy                dofz               dofrx               dofry               dofrz
                   2                                       1                   1                   1                   1                   1                   1
$ load `tip_load`
*LOAD_NODE_SET
$               nsid                 dof                lcid                  sf
                   3                   3                   1               -25.0
*DATABASE_BINARY_D3PLOT
$                 dt
                 1.0
*DATABASE_SPCFORC
$                 dt
                 1.0
*DATABASE_EXTENT_BINARY
$              neiph               neips              maxint              strflg
                                                                               1
*END
