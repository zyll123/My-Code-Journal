import subprocess
import datetime

def sync_to_github():
    # 获取当前精确时间作为 Git 提交的信息
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    commit_msg = f"Auto-sync: {now} 日志更新"
    
    print("🚀 开始自动扫描更改并同步到 GitHub...")
    
    try:
        # 1. 将所有新写的内容或修改添加到暂存区
        subprocess.run(["git", "add", "."], check=True)
        
        # 2. 执行 commit 提交
        # 这一步需要捕获输出，因为如果日记没有改动，git commit 会返回提示
        result = subprocess.run(["git", "commit", "-m", commit_msg], capture_output=True, text=True)
        
        if "nothing to commit" in result.stdout or "无文件要提交" in result.stdout:
            print("✨ 当前日记没有任何修改，无需同步。")
            return
            
        # 3. 推送到 GitHub (假设你的主分支叫 main)
        print("⏳ 正在推送到云端，请稍候...")
        subprocess.run(["git", "push", "origin", "main"], check=True)
        
        print(f"☁️ 同步完美结束！今日的数据已于 {now} 安全上云。")
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ 同步过程中出现错误，请检查网络或 Git 配置。\n错误信息: {e}")

if __name__ == "__main__":
    sync_to_github()