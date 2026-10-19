*KEYWORD long=y
$ Written by FieldCraft: implicit linear analysis, long format (20-column fields).
$ Units: m-kg-s-N-Pa (metre, kilogram, second, newton, pascal, kelvin, consistent).
*CONTROL_UNITS
$             length                time                mass                temp
                   m                 sec                  kg                   k
*NODE
$                nid                   x                   y                   z
                   1                 0.0                 0.0                 0.0
                   2                 1.0                 0.0                 0.0
                   3                 0.0                 1.0                 0.0
                   4                 0.0                 0.0                 1.0
                   5                 2.0                 0.0                 0.0
                   6                 3.0                 0.0                 0.0
                   7                 2.0                 1.0                 0.0
                   8                 2.0                 0.0                 1.0
                   9                 2.5                 0.0                 0.0
                  10                 2.5                 0.5                 0.0
                  11                 2.0                 0.5                 0.0
                  12                 2.0                 0.0                 0.5
                  13                 2.5                 0.0                 0.5
                  14                 2.0                 0.5                 0.5
                  15                 4.0                 0.0                 0.0
                  16                 5.0                 0.0                 0.0
                  17                 4.0                 1.0                 0.0
                  18                 4.0                 0.0                 1.0
                  19                 5.0                 0.0                 1.0
                  20                 4.0                 1.0                 1.0
                  21                 6.0                 0.0                 0.0
                  22                 7.0                 0.0                 0.0
                  23                 7.0                 1.0                 0.0
                  24                 6.0                 1.0                 0.0
                  25                 6.0                 0.0                 1.0
                  26                 7.0                 0.0                 1.0
                  27                 7.0                 1.0                 1.0
                  28                 6.0                 1.0                 1.0
                  29                 6.5                 0.0                 0.0
                  30                 7.0                 0.5                 0.0
                  31                 6.5                 1.0                 0.0
                  32                 6.0                 0.5                 0.0
                  33                 6.5                 0.0                 1.0
                  34                 7.0                 0.5                 1.0
                  35                 6.5                 1.0                 1.0
                  36                 6.0                 0.5                 1.0
                  37                 6.0                 0.0                 0.5
                  38                 7.0                 0.0                 0.5
                  39                 7.0                 1.0                 0.5
                  40                 6.0                 1.0                 0.5
            99999999                 6.5                 0.5                 2.0
*MAT_ELASTIC_TITLE
Aluminium
$                mid                  ro                   e                  pr
                   1              2700.0       70000000000.0                0.33
*SECTION_SOLID
$              secid              elform
                   1                  23
*SECTION_SOLID
$              secid              elform
                   2                  16
*SECTION_SOLID
$              secid              elform
                   3                  13
*SECTION_SOLID
$              secid              elform
                   4                   2
*PART
hex20
$                pid               secid                 mid               eosid                hgid
                   1                   1                   1
*PART
tet10
$                pid               secid                 mid               eosid                hgid
                   2                   2                   1
*PART
tet4
$                pid               secid                 mid               eosid                hgid
                   3                   3                   1
*PART
wedge
$                pid               secid                 mid               eosid                hgid
                   4                   4                   1
*ELEMENT_SOLID
$                eid                 pid
$                 n1                  n2                  n3                  n4                  n5                  n6                  n7                  n8                  n9                 n10
                   1                   3
                   1                   2                   3                   4                   4                   4                   4                   4
                   2                   2
                   5                   6                   7                   8                   9                  10                  11                  13                  14                  12
                   3                   4
                  16                  15                  18                  19                  17                  17                  20                  20
*ELEMENT_SOLID_H20
$                eid                 pid
$                 n1                  n2                  n3                  n4                  n5                  n6                  n7                  n8                  n9                 n10
                   4                   1
                  21                  22                  23                  24                  25                  26                  27                  28                  29                  30
                  31                  32                  33                  34                  35                  36                  37                  38                  39                  40
*SET_NODE_LIST_TITLE
base
$                sid
                   1
                   1                   2                   3                   4
*SET_NODE_LIST_TITLE
lid
$                sid
                   3
                  25                  26                  27                  28                  33                  34                  35                  36
*SET_NODE_LIST_TITLE
tet10
$                sid
                   4
                   5                   6                   7                   8                   9                  10                  11                  12
                  13                  14
$ coupling `spider`
*CONSTRAINED_INTERPOLATION
$               icid                dnid                ddof                cidd                ityp               idnsw                 fgm
                   1            99999999              123456                                       1                                       1
$               inid                idof
                   3                 123
*TITLE
Load
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
Load ramp
$               lcid
                   1
$                 a1                  o1
                 0.0                 0.0
                 1.0                 1.0
*BOUNDARY_SPC_SET_ID
                   1hold
$               nsid                 cid                dofx                dofy                dofz               dofrx               dofry               dofrz
                   1                                       1                   1                   1                   1                   1                   1
*BOUNDARY_PRESCRIBED_MOTION_SET_ID
                   2lift
$               nsid                 dof                 vad                lcid                  sf
                   4                   3                   2                   1               0.001
*LOAD_NODE_POINT
$                nid                 dof                lcid                  sf
                  25                   1                   1  -83.33333333333337
                  26                   1                   1  -83.33333333333336
                  27                   1                   1  -83.33333333333323
                  28                   1                   1   -83.3333333333333
                  33                   1                   1   333.3333333333333
                  34                   1                   1  333.33333333333337
                  35                   1                   1   333.3333333333334
                  36                   1                   1   333.3333333333333
            99999999                   3                   1                50.0
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
