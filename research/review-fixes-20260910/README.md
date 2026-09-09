# 代码审查修复记录（2026-09-10）

## 已修复

- dynamic统一重置位置索引，校验重复/缺失个体年份；未配对个体保留缺失结果。缓存观测与参照期距离，不改变定义。
- run模块导入不创建目录；新增configure_output及--output-dir，安装环境默认输出到当前工作目录。
- strict=True遇到不可行异常抛RuntimeError；strict=False返回失败状态，求解器保持可用于后续模型。
- 最低Python版本统一为3.9；wheel仅标记Python3；新增book/test可选依赖。
- 正式案例统一通过loaders.read_raw和deabook.data_path读取五份数据，移除用户缓存自动注入。
- 第13章删除未使用的额外求解，每次复算少30个LP。
- unique_cases.py和8份未用CSV移入research/legacy-case-candidates-20260910，不随安装包发布。

## 验证

- tests/：50 passed、3 skipped；新增17项回归覆盖索引不变性、手算指数、重复/未配对数据、距离缓存、输出位置、数据入口以及求解失败后恢复。
- 3项跳过为CNLS外部求解器和amplpy缺失，不视为通过。
- 13章主案例完整重跑并与原full.csv对比：全部通过，rtol=1e-5、atol=1e-7。正式结果与图表未被验证输出覆盖。
- pip wheel成功；核对Python>=3.9元数据和五份正式数据；在临时目录解包wheel，使用wheel中的代码读取全部数据并通过CLI生成第1章图表，验证包内代码不依赖书稿目录或源码路径。第三方依赖使用当前验证环境，未进行Python3.9独立环境安装测试。

证据：tests.txt、case-rerun.txt、build.txt、wheel-check.txt。复算脚本：run_cases.py、check_wheel.py。before/为本轮修改前备份，edit.py仅为一次性历史编辑记录，不应重跑。

本次仅修改和验证本地文件，未推送GitHub。

## MOSEK补充验证

用户确认本机已有许可证后，检查发现默认许可证文件存在，原验证环境仅缺MOSEK Python接口。补装接口后使用MOSEK 11.2.3实际求解CNLS，终止状态optimal，12个残差均为有限值。全套测试重新运行：51 passed、2 skipped；剩余两项仅因amplpy缺失。证据：tests-with-mosek.txt、mosek-check.txt；复现脚本check_mosek.py。未读取或复制许可证内容。此前“缺少CNLS求解器”的结论仅适用于补装前的验证环境，不能解释为本机没有许可证。
