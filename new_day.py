import os
import datetime

# 1. 获取今天的日期，作为文件名（例如：2026-05-15.md）
today = datetime.date.today().strftime("%Y-%m-%d")
file_name = f"{today}.md"

# 使用变量定义 Markdown 的代码块符号，防止渲染冲突导致代码无法复制
md_code = "```"

# 2. 定义你的专属日记模板
template = f"""# {today} 研发与创作者日志

## 🧠 Creator's Flow (灵感与架构)
- 

## 🔬 Research & Code (生信分析 / 代码沉淀)
### 核心代码片段备份
{md_code}python
# 在这里记录有用的 Python / R (Seurat/Scanpy) 代码片段

{md_code}

### 分析复盘与文献阅读
- 

## 💼 Career & Life (职业规划与随笔)
- 
"""

# 3. 创建并写入当天的 Markdown 文件
if not os.path.exists(file_name):
    # 使用 utf-8 编码写入，防止中文字符乱码
    with open(file_name, "w", encoding="utf-8") as file:
        file.write(template)
    print(f"✅ 成功生成今日日记模板: {file_name}")
else:
    print(f"⚠️ 文件 {file_name} 已存在，请直接在 VS Code 中继续记录！")

# 4. 自动在 VS Code 中打开这个生成好的文件
os.system(f"code {file_name}")