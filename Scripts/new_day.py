import os
import datetime

# ==========================================
# ⚙️ 配置区：在这里修改你的路径
# ==========================================
# 相对路径写法：因为脚本在 Scripts 里，"../" 表示退回上一级目录
TARGET_DIR = "../00_Daily_Drafts" 

# 如果你想用绝对路径（直接指定死位置），可以改成类似下面这样：
# TARGET_DIR = "D:/AI/My-Code-Journal/00_Daily_Drafts"
# ==========================================

# 1. 获取今天的日期，并生成最终的文件保存路径
today = datetime.date.today().strftime("%Y-%m-%d")
file_name = f"{today}.md"

# 🌟 关键点 1：确保目标文件夹存在，如果不存在，Python 会自动帮你建一个
os.makedirs(TARGET_DIR, exist_ok=True)

# 🌟 关键点 2：将文件夹路径和文件名安全地拼接在一起
full_path = os.path.join(TARGET_DIR, file_name)

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

# 3. 创建并写入当天的 Markdown 文件（注意这里统一换成了 full_path）
if not os.path.exists(full_path):
    with open(full_path, "w", encoding="utf-8") as file:
        file.write(template)
    print(f"✅ 成功生成今日日记模板: {full_path}")
else:
    print(f"⚠️ 文件 {file_name} 已存在，请直接在 VS Code 中继续记录！")

# 4. 自动在 VS Code 中打开这个生成好的文件
os.system(f"code {full_path}")