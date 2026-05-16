

# --- Merged via AI on 2026-05-16 ---
# Differential GRN Analysis: AD vs WT
---
**Goal:** For each cell type × month, identify regulons with significantly different activity between AD and WT (Mann-Whitney U test), filtered to regulons with ≥10% activation rate in at least one group.

**Input:**
- `combined`: Seurat object with metadata (genotype, month, class_id_label)
- `activity`: Binarized TF activity matrix (0/1 per cell)
- `AUC_gene`: Regulon AUC activity matrix (continuous scores)

# --- Merged via AI on 2026-05-16 ---
## 4. Differential Analysis
For each **cell type × month** combination:
1. Filter regulons ≥10% activation in AD **or** WT
2. Mann-Whitney U test on AUC scores
3. BH-FDR correction

# --- Merged via AI on 2026-05-16 ---
## 6. Cross-Brain-Region Comparison
Only runs when `split_by_brain <- TRUE`. Parses the `CellType` column to extract base cell type and brain region, then compares which regulons are shared vs brain-region-specific.

# --- Merged via AI on 2026-05-16 ---
🌟 阿兹海默症的"破局者"！深理工叶克强教授全科普 🧠
今天必须安利一位神经科学领域的"真·大佬"——叶克强教授！他不仅发文章发到手软，更是真正把科研变成救命药的人。看完你就知道什么叫做科学家该有的样子 👇
✦ ✦ ✦
🌍 学术履历有多硬核？
20 年埃默里大学（Emory）执教经验，却毅然全职回国。现在是深圳理工大学生命健康学院讲席教授 + 中科院深圳先进院研究员。几十年只做一件事：攻破阿尔茨海默病、帕金森病等神经退行性疾病。学术天花板，就是这么顶 🏆
✦ ✦ ✦
🔬 两大颠覆性发现，每一个都够发 Cell
1️⃣ 锁定致病"元凶" —— Tau N368 片段
Tau 蛋白是 AD 研究的关键，叶教授团队发现天冬酰胺内肽酶（AEP）剪切产生的 Tau N368 片段，既是早期诊断的金标志物，更是药物靶点的"靶心" 🎯
2️⃣ 肠脑轴的奇妙证据 🦠➡️🧠
团队在《Cell》子刊证实：Tau 和 α-突触核蛋白的病理竟能从肠道"顺藤摸瓜"传到大脑！从此神经退行性疾病多了一个全新视角
✦ ✦ ✦
📚 代表性文献速览
▎Nature Medicine（2014）
⭐ 发现 δ-secretase 剪切 Tau N368 → 神经缠结
▎Nature（2022）
阻断 FSH 改善 AD 小鼠认知，解释女性高发机制
▎Nat Struct Mol Biol（2022）
DOPEGAL 修饰 Tau 促聚集，揭示蓝斑核易损之谜
▎Cell（2023）
全球首个 α-突触核蛋白 PET 示踪剂
▎Nat Commun ×4（2015–23）
δ-secretase 抑制剂、C/EBPβ 通路、FSH+ApoE4
📌 核心学术主线：C/EBPβ → δ-secretase → Tau N368，一条从机制到药物的完整链条
✦ ✦ ✦
💊 不只是发论文，他做出了真药
教授 + 创业者双重身份，2015 年创办 博芮健制药，落户深圳脑创中心。核心管线 BrAD-R13 片 是全球首款进入临床的口服 TrkB 受体激动剂，精准激活脑内保护通路，从源头减少致病蛋白。I 期临床已顺利收官，给全球千万 AD 患者带来真实可触的希望 🌅
✦ ✦ ✦
👨‍🔬 为什么他值得被更多人知道？
做基础研究的人很多，但敢跨越"死亡之谷"、把论文变成药片的人太少。从上海到深圳，叶教授一路追寻最好的转化生态，用硬核成果证明了中国科学家的担当 💪
✦ ✦ ✦
让我们一起期待创新药早日上市，造福千万患者 ❤️
✦ ✦ ✦
#阿尔茨海默症 #神经科学 #科研日常 #药学生物 #深圳理工大学 #中科院深圳先进院 #创新药 #Tau蛋白 #脑科学 #科学家故事