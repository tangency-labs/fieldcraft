*KEYWORD long=y
$ Written by FieldCraft: implicit linear analysis, long format (20-column fields).
$ Units: mm-t-s-N-MPa (millimetre, tonne, second, newton, megapascal, degree Celsius, consistent).
*CONTROL_UNITS
$             length                time                mass                temp
                  mm                 sec            mtrc_ton                   c
*NODE
$                nid                   x                   y                   z
                   1                 0.0                 0.0                 0.0
                   2              1000.0                 0.0                 0.0
                   3              1000.0                 0.0               500.0
                   4                 0.0                 0.0               500.0
*MAT_ELASTIC_TITLE
Steel
$                mid                  ro                   e                  pr
                   1             7.85e-9            210000.0                 0.3
*SECTION_BEAM
$              secid              elform
                   1                   3
$                  a
                50.0
*SECTION_BEAM
$              secid              elform
                   2                   2
$                  a                 iss                 itt                   j                  sa
              1200.0            360000.0             40000.0  126392.12692253345  1019.6078431372549
*SECTION_BEAM
$              secid              elform
                   3                   2
$                  a                 iss                 itt                   j                  sa
               444.0            101972.0            101972.0            151959.0               222.0
*SECTION_BEAM
$              secid              elform
                   4                   2
$                  a                 iss                 itt                   j                  sa
   314.1592653589793   7853.981633974483   7853.981633974483  15707.963267948966   278.4593488409135
*SECTION_BEAM
$              secid              elform
                   5                   2
$                  a                 iss                 itt                   j                  sa
  238.76104167282426   43215.74854278119   43215.74854278119   86431.49708556238   127.2523276188772
*PART
brace
$                pid               secid                 mid               eosid                hgid
                   1                   1                   1
*PART
girder
$                pid               secid                 mid               eosid                hgid
                   2                   2                   1
*PART
leg
$                pid               secid                 mid               eosid                hgid
                   3                   3                   1
*PART
post
$                pid               secid                 mid               eosid                hgid
                   4                   4                   1
*PART
rail
$                pid               secid                 mid               eosid                hgid
                   5                   5                   1
*ELEMENT_BEAM
$                eid                 pid                  n1                  n2                  n3                  n4                  n5                  n6                  n7                  n8
                   5                   1                   1                   3
*ELEMENT_BEAM_ORIENTATION
$                eid                 pid                  n1                  n2                  n3
$                 vx                  vy                  vz
                   1                   2                   1                   2
                 0.0                 1.0                 1.0
                   2                   4                   2                   3
                 1.0                 0.0                 0.0
                   3                   5                   3                   4
                 0.0                 0.0                 1.0
                   4                   3                   4                   1
                 1.0                 0.0                 0.0
*SET_NODE_LIST_TITLE
base
$                sid
                   1
                   1
*SET_NODE_LIST_TITLE
corner
$                sid
                   3
                   3
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
                   1clamp
$               nsid                 cid                dofx                dofy                dofz               dofrx               dofry               dofrz
                   1                                       1                   1                   1                   1                   1                   1
$ load `push`
*LOAD_NODE_SET
$               nsid                 dof                lcid                  sf
                   3                   1                   1               500.0
                   3                   3                   1             -1000.0
$ load `turn`
*LOAD_NODE_SET
$               nsid                 dof                lcid                  sf
                   3                   6                   1            200000.0
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
