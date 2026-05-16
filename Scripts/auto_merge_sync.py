import os
import re
import datetime
import subprocess
from openai import OpenAI

# ==========================================
# ⚙️ 第一部分：AI 路由器配置
# ==========================================

# 1. 填入你的大模型 API 密钥 (这里以 DeepSeek 为例，因为它便宜且国内直连)
# 只要是支持 OpenAI SDK 的模型都可以无缝替换 (如 Kimi, Qwen, ChatGPT)
API_KEY = "sk-b0e56f5da20648fd99b8fabbd09bd7d9" 
BASE_URL = "https://api.deepseek.com" # 如果用其他模型，替换为对应的官方 URL

client = OpenAI(api_key=API_KEY, base_url=BASE_URL)

# 2. 定义你的“知识图谱架构” (AI 将根据这些键名进行分类)
KNOWLEDGE_BASE = {
    "Seurat_R_Pipeline": "../01_Bioinformatics/Seurat_Recipes.R",
    "Scanpy_Python_Pipeline": "../01_Bioinformatics/Scanpy_Pipelines.py",
    "AD_Pathology_Notes": "../02_Research_Notes/AD_Analysis_Tools.md",
    "Career_MNC_Interview": "../03_Career_Track/Interview_Notes.md"
}

DRAFTS_DIR = "../00_Daily_Drafts"

# ==========================================
# 🧠 第二部分：AI 核心分析逻辑
# ==========================================

def ask_ai_to_classify(content_snippet):
    """让大模型阅读代码/文本，并决定它属于哪个分类"""
    
    categories_str = ", ".join(KNOWLEDGE_BASE.keys())
    
    system_prompt = f"""
    你是一个顶级的医学与生信科研助手。你的任务是对用户提供的代码或研究笔记进行精确分类。
    请阅读用户输入的内容，并将其归类到以下【唯一个】类别中：
    [{categories_str}]
    
    判断标准提示：
    - 如果是 R 语言处理单细胞数据，通常属于 Seurat_R_Pipeline
    - 如果是 Python 处理单细胞或深度学习，通常属于 Scanpy_Python_Pipeline
    - 如果是关于阿尔茨海默症(AD)、微胶质细胞、靶点网络，属于 AD_Pathology_Notes
    - 如果是关于外企(MNC)面试、药企职位、简历，属于 Career_MNC_Interview
    
    如果内容非常混乱或不属于任何分类，请输出 "Unknown"。
    【绝对规则】：你的输出只能是上述类别名称或 "Unknown"，不要有任何多余的标点、换行或解释！
    """

    try:
        response = client.chat.completions.create(
            model="deepseek-chat", # 替换为你使用的具体模型名称
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": content_snippet}
            ],
            temperature=0.1, # 极低的温度保证分类的确定性和稳定性
            max_tokens=10
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        print(f"⚠️ AI 分类服务请求失败: {e}")
        return "Unknown"

# ==========================================
# 🚀 第三部分：提取、分发与同步逻辑
# ==========================================

# ==========================================
# 🚀 第三部分：提取、分发与同步逻辑
# ==========================================

def extract_and_route(file_path):
    """读取草稿，提取所有内容块并交给 AI 分发 (带自动去重)"""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # 提取 Markdown 中的所有代码块
    code_blocks = re.findall(r'```([a-zA-Z]*)\n(.*?)```', content, re.DOTALL)
    
    merged_count = 0
    for lang, code_text in code_blocks:
        code_text = code_text.strip()
        
        # 【新增防御 1】跳过空代码，并过滤掉日记模板自带的占位提示词
        if not code_text or "在这里记录有用的" in code_text:
            continue
            
        print(f"🤖 AI 正在阅读一段 {lang or '未知语言'} 代码...")
        
        # 将代码喂给 AI 判断
        target_category = ask_ai_to_classify(f"语言:{lang}\n内容:\n{code_text}")
        
        if target_category in KNOWLEDGE_BASE:
            target_file = KNOWLEDGE_BASE[target_category]
            
            # 【新增防御 2】查重机制：判断这段代码是否已经存在于目标文件中
            is_duplicate = False
            if os.path.exists(target_file):
                with open(target_file, 'r', encoding='utf-8') as check_f:
                    existing_content = check_f.read()
                    if code_text in existing_content:
                        is_duplicate = True
            
            if is_duplicate:
                print(f"⏩ 代码查重: 发现完全相同的代码，已自动跳过合并。")
                continue  # 如果重复，直接跳过，进入下一段代码
            
            # 如果不重复，则正常追加写入
            os.makedirs(os.path.dirname(target_file), exist_ok=True)
            with open(target_file, 'a', encoding='utf-8') as target_f:
                target_f.write(f"\n\n# --- Merged via AI on {datetime.date.today()} ---\n")
                target_f.write(code_text)
            
            print(f"🎯 AI 判定归类: 【{target_category}】 -> 已合并至 {os.path.basename(target_file)}")
            merged_count += 1
        else:
            print(f"🤷 AI 认为这段内容属于: 【{target_category}】(丢弃或手动处理)")
                
    return merged_count

def sync_to_github():
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    print("\n🚀 开始自动扫描更改并同步到 GitHub...")
    try:
        os.chdir("..") 
        subprocess.run(["git", "add", "."], check=True)
        result = subprocess.run(["git", "commit", "-m", f"AI-Merge & Sync: {now}"], capture_output=True, text=True)
        
        if "nothing to commit" in result.stdout or "无文件要提交" in result.stdout:
            print("✨ 没有任何修改，无需同步。")
        else:
            print("⏳ 正在通过 SSH 通道推送到云端...")
            subprocess.run(["git", "push", "origin", "main"], check=True)
            print(f"☁️ 完美收工！你的数据图谱已自动构建并于 {now} 备份。")
    except Exception as e:
        print(f"\n❌ 同步失败: {e}")

if __name__ == "__main__":
    today_draft = os.path.join(DRAFTS_DIR, f"{datetime.date.today()}.md")
    
    if os.path.exists(today_draft):
        print(f"🔍 正在扫描今日草稿: {os.path.basename(today_draft)}")
        count = extract_and_route(today_draft)
        if count == 0:
            print("📝 今日草稿中没有识别出有效的代码块。")
    else:
        print("📭 今天没有草稿。")
        
    sync_to_github()