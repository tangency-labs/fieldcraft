*KEYWORD long=y
$ Written by FieldCraft: implicit linear analysis, long format (20-column fields).
$ Units: mm-t-s-N-MPa (millimetre, tonne, second, newton, megapascal, degree Celsius, consistent).
*CONTROL_UNITS
$             length                time                mass                temp
                  mm                 sec            mtrc_ton                   c
*NODE
$                nid                   x                   y                   z
                   1                 0.0                 0.0                 0.0
                   2                10.0                 0.0                 0.0
                   3                20.0                 0.0                 0.0
                   4                 0.0                10.0                 0.0
                   5                10.0                10.0                 0.0
                   6                20.0                10.0                 0.0
                   7                30.0                 0.0                 0.0
                   8                30.0                10.0                 0.0
*MAT_ELASTIC_TITLE
Steel
$                mid                  ro                   e                  pr
                   1             7.85e-9            210000.0                 0.3
*SECTION_SHELL
$              secid              elform                shrf                 nip
                   1                   4  0.8333333333333334                   5
$                 t1                  t2                  t3                  t4                nloc
                 0.5                                                                             0.5
*SECTION_SHELL
$              secid              elform                shrf                 nip
                   2                  16  0.8333333333333334                   5
$                 t1                  t2                  t3                  t4                nloc
                 0.5                                                                             0.5
*PART
skin tri3
$                pid               secid                 mid               eosid                hgid
                   1                   1                   1
*PART
skin quad4
$                pid               secid                 mid               eosid                hgid
                   2                   2                   1
*ELEMENT_SHELL
$                eid                 pid                  n1                  n2                  n3                  n4                  n5                  n6                  n7                  n8
                   1                   2                   1                   2                   5                   4
                   2                   2                   2                   3                   6                   5
                   3                   1                   3                   7                   8                   8
                   4                   1                   3                   8                   6                   6
*SET_NODE_LIST_TITLE
edge
$                sid
                   1
                   1                   4
*SET_SEGMENT_TITLE
top
$                sid
                   3
$                 n1                  n2                  n3                  n4
                   1                   2                   5                   4
                   2                   3                   6                   5
                   3                   7                   8                   8
                   3                   8                   6                   6
*TITLE
Sag
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
Sag ramp
$               lcid
                   1
$                 a1                  o1
                 0.0                 0.0
                 1.0                 1.0
*BOUNDARY_SPC_SET_ID
                   1hold
$               nsid                 cid                dofx                dofy                dofz               dofrx               dofry               dofrz
                   1                                       1                   1                   1                   1                   1                   1
$ load `press`
*LOAD_SEGMENT_SET
$               ssid                lcid                  sf
                   3                   1                0.01
*LOAD_BODY_Z
$               lcid                  sf
                   1              9810.0
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
