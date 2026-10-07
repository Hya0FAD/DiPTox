# 条件列单位与对数转换 / Condition column transformations

可用于暴露时长、温度、浓度等带单位的条件列。条件列转换和单组合并自 1.1.2 提供；JSON 数组批量合并自 1.1.3 提供。API、GUI、CLI 均支持。

## GUI

在单位转换页面将 **Apply to** 设为 **Condition column**，选择 **Value column** 和 **Unit column**。
填写 **Standard unit**，并在 **Condition transformation** 选择 `None`、`log10` 或 `-log10`。
换算先执行，对数以 10 为底。只做对数时清空 Standard unit；只换单位时选择 None。
自定义公式及分子量来源沿用同页设置。完成后在去重页选择生成的条件列。
**Filter by column values** 和 **Merge column values** 位于独立的 **Column adjustments（列调整）** 页面，侧栏顺序在去重上方，可先整理参数，再分组去重。

## 合并多个类别为一个值

在列调整页的 **Merge column values** 卡片选择列，再多选待合并值（每个值附带行数），
在 **Replace with** 输入 `other` 等标签。可自定义输出列名，默认是 `原列名 (Merged)`。
未选择的值保持原样，所有行保留且不会进入排除表，原始列也保留。完成后在去重条件中选择新列。
空选择不做修改，缺失值可单独选中合并；支持撤销。同名输出列会报错，可换新列名或撤销上一次操作。
要一次完成多组分类，在 GUI 选择 **Multiple groups (JSON)**，粘贴规则数组后点击 **Apply merge**。
只生成一个输出列，所有规则同时匹配原始值，不会级联替换。未列出的值保留原样；同一原值对应不同标签会报错，且不修改数据。
空数组 `[]` 不做修改，每组 `values` 必须非空。匹配区分大小写和空格。
完整生命阶段分类数组见 [life_stage_merge_groups.json](../examples/life_stage_merge_groups.json)。

```json
[
  {"values": ["Embryo", "Egg", "Blastula"], "replacement": "Embryonic"},
  {"values": ["Larva", "Fry", "Alevin"], "replacement": "Larval / post hatch"},
  {"values": ["Sperm", "Oocyte / ova"], "replacement": "Gamete"}
]
```

API 使用 `pipeline.merge_column_values('Stage', groups=rules, output_column='Stage group')`。
CLI 使用 `merge-values --column Stage --groups '<JSON数组>'`，同时提供 `--input`、`--output`。
流水线配置把数组放在 `merge-values` 步骤的 `params.groups`；`groups` 与单组 `values` 不可同时指定。
原有单组功能仍可使用。

```python
group = pipeline.merge_column_values(
    'Species', ['mouse', 'rabbit'], replacement='other',
    output_column='Species group',
)
pipeline.config_deduplicator(condition_cols=[group])
pipeline.dataset_deduplicate()
```

CLI 的 `merge-values` 接受 `--column`、`--values`（JSON 数组）、`--replacement` 和可选 `--output-column`。
JSON 流水线示例：

```json
{"op": "merge-values", "params": {
  "column": "Species", "values": ["mouse", "rabbit"],
  "replacement": "other", "output_column": "Species group"
}}
```

API/JSON 中数字 `1`、文本 `"1"`、布尔值 `true` 独立匹配；`null` 匹配缺失单元格。
替换标签为非空文本。完整示例见 `UseExample.py` 的 `use_case_8()`。

## API

```python
value_col, unit_col = pipeline.transform_column(
    value_col='Duration', unit_col='Duration unit',
    standard_unit='h', log_transform='-log10',
)
pipeline.config_deduplicator(condition_cols=[value_col, unit_col])
pipeline.dataset_deduplicate()
```

`1 d` 和 `24 h` 先统一为 `24 h`，再转换为约 `-1.380211`，单位标记为 `-log10(h)`。
仅转换单位的输出列以 ` (Standardized)` 结尾；包含对数的以 ` (Transformed)` 结尾。
返回值是生成的数值列和单位列名。原始列及主要目标列映射保留，推荐导出包含生成的数值/单位列。

`conversion_rules={('week', 'h'): 'x * 168'}` 可补充换算规则。需要分子量时显式指定
`molecular_weight_source='original'/'standardized'`，或 `molecular_weight_col='MW'`；标准化来源要求先预处理。
可以逐次处理多对列。输出重名时请撤销上次操作，或选择另一来源列。

## CLI

```shell
diptox transform-column --input input.csv --value-col Duration --column-unit "Duration unit" --standard-unit h --log-transform=-log10 --output result.csv --excluded excluded.csv
```

`--column-unit` 指本次条件单位列；原有 `--unit-col` 仍用于映射主要目标的单位列。
负对数参数建议使用 `--log-transform=-log10` 写法。自定义换算通过 `--conversion-rules rules.json` 传入。

流水线 JSON 中对应步骤为：

```json
{"op": "transform-column", "params": {
  "value_col": "Duration", "unit_col": "Duration unit",
  "standard_unit": "h", "log_transform": "-log10"
}}
```

## 数据处理规则

- 自 1.1.3 起，单位换算不按原始有效数字舍入，后续对数与聚合使用完整浮点结果。需要显示少量小数位时应使用显示格式，不要把舍入值写回计算列。旧的已舍入结果需从原值重新计算。
- 缺失/非数值/非有限数值、无法换算的行进入排除表并记录原因；取对数还要求数值大于零。
- 缺少整套换算规则、列不存在或参数错误会终止操作，API 和 GUI 不提交部分结果。
- 已标记为对数尺度的列不能再次进行线性换算或取对数，请选择原始线性数值。
- 仅取对数不统一单位；若存在多个单位，去重条件必须同时包含输出数值列和单位列，或先统一单位。
- `pipeline.undo()` 撤销转换及其排除记录。负对数单位通过现有 Excel 安全文本导出逻辑保存。

## English

Select **Condition column** on the unit page and choose a value/unit pair. Unit conversion runs first,
followed by optional base-10 logarithms. Leave the standard unit blank for log-only processing.
The API returns generated value/unit column names; original data and primary target mappings remain intact.
Use both returned columns as deduplication conditions, especially when log-only processing leaves mixed units.
Invalid rows are excluded with reasons, and undo restores them. Explicit molecular-weight selection is required
for mass/molar conversions. See `use_case_7()` in `UseExample.py` for a complete offline example.
