# Route 8 面支承接口

FieldCraft 面支承命令：

```json
{"command":"model.foundation.create","params":{"name":"BASE_Bed","set":"bottom","normal":6e9,"tangential":2e9,"normal_damping":0,"tangential_damping":0}}
```

系数单位为 N/m³；阻尼为 N·s/m³。10 mm BASE（91449 节点）用于两点静力 E+k 识别：顶面中心 `(-163,2,21.9) mm`、正下方底面 `(-163,2,-30) mm`，均取竖向位移，真值 E=230 GPa、kn=6e9 N/m³、kt=2e9 N/m³。协调者记录的识别为442 s、10次评估、误差6.5e-9/1.1e-9、秩2、条件数686。

5 mm 仅用于一次面支承模态频率证明，不进入识别循环；前10阶为505.38、513.70、531.84、880.05、917.45、959.86、1091.22、1232.41、1746.45、1886.06 Hz，均在0–2000 Hz。5 mm单次静力成本221 s。该证据来自本路记录，待项目验收；不代表商业软件结果。
