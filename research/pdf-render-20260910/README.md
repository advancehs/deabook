# 整书 PDF 渲染记录

2026-09-10：使用 RStudio 随附 Quarto 1.8.25 与 TeX Live 2026 LuaLaTeX 完成整书 PDF 渲染，退出码 0。

- 输出：docs/deabook文档/_book/资源环境绩效评估方法：DEA模型与Python应用.pdf
- 同步更新书稿目录下的同名 PDF；旧 PDF 保存在 previous-book.pdf。
- 最终 464 页，4,412,255 字节。
- 修正默认字体缺字：正文 Times New Roman，代码 Consolas，星形符号 Segoe UI Symbol。这些字体需在渲染机器上安装。
- Pandoc/Quarto 三轮 LuaLaTeX 编译完成。
- 用 PyMuPDF 与 pdftotext 核对最终 PDF，定位第12章正文313、314、325页公式。
- 未完成整书逐页视觉审校；首次编译发现其他长公式及部分内容超宽，不能将编译通过视为全部排版合格。
- render.py 为本机 Windows Python 启动脚本，补足从 WSL 启动时的 Windows 环境变量。
