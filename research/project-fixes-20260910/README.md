# 全项目审查修复（2026-09-10）

本轮修复了项目审查中的七项问题，并补充同类DDF边界及网站示例修复。

- 径向theta/phi允许合法跨期范围；DDF允许负跨期方向距离，手算MQ=2与径向一致。
- 核心数据入口拒绝重复DataFrame索引；面板入口验证非空、非缺失id/year、唯一索引和唯一id-year。
- CNLS optimize_model转用core严格求解状态检查，移除裸except自动回退；MOSEK许可证探测与测试失败区分。
- 跨期距离框保留时期元数据，分组shift只接受面板相邻期；BPG派生框同步元数据，缺期结果保留NaN。
- 正整数窗口严格验证；偶数窗口向未来多分配一期，边界截断。
- 网站Notebook使用CONTEMPORARY，修复重复面板构造及pivot调用；测试检查最终6×2有限值结果。
- Linux/Windows/macOS CI安装.[test]并运行pytest，移除GDAL；文档CI改pytest，PyPI发布前增加测试。
- mqrd分解改为VRS跨期技术变化与单独SCALECH_RD，满足MQ=EFFCH*TECHCH_RD*SCALECH_RD；更新第4章公式和接口说明。验证依据是距离函数分解与手算代数测试，本轮未重新取得原期刊全文，不能称为原文逐项复核。

## 验证

71 passed、2 skipped（仅amplpy缺失）；含新增20项回归测试。MOSEK CNLS测试通过。
13章主案例全部重跑并与原完整结果比较通过（rtol=1e-5,atol=1e-7）；正式图表未被覆盖。
wheel重新构建并在临时目录验证五份数据读取和第1章CLI图表生成。
工作流YAML静态检查；GitHub Actions及整站MkDocs尚未在远端实际执行。网站Notebook代码已在本地顺序执行测试。

证据：tests.txt、cases.txt、build.txt、wheel.txt；tests/test_project_regressions.py为新增回归。edit.py仅记录一次性历史编辑，不要重跑。

未推送GitHub。此次跨期和mqrd包接口结果有意修正，旧包接口生成的动态表应重新计算；当前book_cases主案例结果未变。
