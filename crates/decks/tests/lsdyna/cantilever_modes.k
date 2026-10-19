*KEYWORD long=y
$ Written by FieldCraft: implicit linear analysis, long format (20-column fields).
$ Units: mm-t-s-N-MPa (millimetre, tonne, second, newton, megapascal, degree Celsius, consistent).
*CONTROL_UNITS
$             length                time                mass                temp
                  mm                 sec            mtrc_ton                   c
$ 2 steps are written as *CASE analyses: run LS-DYNA with CASE on its command line.
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
            99999999                 5.0                 0.5                 0.5
*MAT_ORTHOTROPIC_ELASTIC_TITLE
Ply
$                mid                  ro                  ea                  eb                  ec                prba                prca                prcb
                   1              1.6e-9            140000.0             10000.0             10000.0 0.02142857142857143 0.02142857142857143                0.45
$                gab                 gbc                 gca                aopt
              5000.0              3500.0              5000.0                 2.0
$                 xp                  yp                  zp                  a1                  a2                  a3
                                                                             1.0                 0.0                 0.0
$                 v1                  v2                  v3                  d1                  d2                  d3
                                                                             0.0                 1.0                 0.0
*SECTION_SOLID
$              secid              elform
                   1                 -18
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
pad (rigid body)
$                sid
                   4
            99999999                   5                  10                  15                  20
*CONSTRAINED_NODAL_RIGID_BODY_TITLE
pad
$                pid                 cid                nsid
                   2                                       4
*CASE_BEGIN_1
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
*LOAD_NODE_POINT
$                nid                 dof                lcid                  sf
            99999999                   2                   1               -10.0
            99999999                   5                   1               250.0
*DATABASE_BINARY_D3PLOT
$                 dt
                 1.0
*DATABASE_SPCFORC
$                 dt
                 1.0
*DATABASE_EXTENT_BINARY
$              neiph               neips              maxint              strflg
                                                                               1
*CASE_END_1
*CASE_BEGIN_2
*TITLE
Modes
*CONTROL_IMPLICIT_GENERAL
$             imflag                 dt0
                   1                 1.0
*CONTROL_IMPLICIT_EIGENVALUE
$               neig              center               lflag              lftend               rflag              rhtend
                   6                                       1                 1.0                   1              5000.0
*CONTROL_TERMINATION
$             endtim
                 1.0
*BOUNDARY_SPC_SET_ID
                   2clamp
$               nsid                 cid                dofx                dofy                dofz               dofrx               dofry               dofrz
                   2                                       1                   1                   1                   1                   1                   1
*CASE_END_2
*END
