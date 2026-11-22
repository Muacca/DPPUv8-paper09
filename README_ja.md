# 固定電荷 Nieh–Yan 半ワームホールの球面喉におけるフェルミオン

Muacca · DPPUv8 Paper 09 · **v1** · 初版 2026-11-22

## 要旨

丸いユークリッド・アクシオン喉では、流束の正負両枝で古典作用が同一でも、異なるフェルミオン状態が実現し得る。空間切片をスピン構造を備えた球面空間形として、Giddings–Strominger 解の反射面から有限の外端までを考える。質量ゼロの試験 Weyl 場では、外端の Atiyah–Patodi–Singer 真空と、最初の準位交差をわずかに越えた側で準備したローレンツ瞬時真空との共分散の差は有限階数である。そのトレースは $`s m_s`$ となる。ここで $`s`$ は Nieh–Yan 電荷の符号、 $`m_s`$ はリフトに関して接触殻に残る Killing スピノルの数である。実定数を係数とする Holst 項を加えると、準位交差が喉から離れ、その間で準備した状態では差は零になる。対応する有限カラー上で、15個のスピンリフトに対する四成分 Dirac 行列式の絶対値の相対対数は $`-2J\eta_Y(0)`$ であり、 $`J`$ は軸性シフトの積分である。無限遠に殻ごとに拡張する際の境界での定義域と極限操作の順序も明示する。実三形式の積分路では、流束を固定したボソン Hessian の負指数は零である。外端の場を固定して連続流束を積分すると、有限半径で斉次かつ偶の負方向が一つ加わり、短いカラーでは Robin 境界項も寄与する。これらは指定した境界実験の条件下における、状態、行列式および安定性の比較である。

[英語版 PDF](DPPUv8-paper09.pdf)・[英語 LaTeX 原稿](LaTeX/main.tex)・[再現資料](script/README.md)・[参考文献](REFERENCES.md)・[履歴](CHANGELOG.md)

## 本文

- [第1節　同じ喉で比べる二つの境界実験](DPPUv8-paper09_sec01.md)
- [第2節　固定電荷と二つの実断面](DPPUv8-paper09_sec02.md)
- [第3節　有限の外端における Morse 指数](DPPUv8-paper09_sec03.md)
- [第4節　喉の境界値とフェルミオンの占有](DPPUv8-paper09_sec04.md)
- [第5節　カラー上の相対 Dirac 行列式の絶対値](DPPUv8-paper09_sec05.md)
- [第6節　境界位相と指定した懸垂](DPPUv8-paper09_sec06.md)
- [第7節　考察](DPPUv8-paper09_sec07.md)

## 付録

- [付録A　規約、流束と向きを指定したレンズ空間の例](DPPUv8-paper09_appA.md)
- [付録B　Robin 条件の符号判定、角モードの比較、スカラー場の積分路](DPPUv8-paper09_appB.md)
- [付録C　Kantowski–Sachs 幾何との比較と曲率による障害](DPPUv8-paper09_appC.md)
- [付録D　コンパクトなゼータ関数の差の角方向解析接続](DPPUv8-paper09_appD.md)
- [付録E　半直線上での閾値の処方](DPPUv8-paper09_appE.md)
- [付録F　比較のための作用と計算の再現](DPPUv8-paper09_appF.md)

## 補足資料の入手と利用

付属リポジトリには、角モードの射影、行列式の絶対値比、Robin 境界条件での指数の検算、および図1の作図ソースを収録した。表の数値を独立に確かめる手順は[付録F](DPPUv8-paper09_appF.md)に記す。[データ案内](data/README.md)・[作図・検算用プログラム](script/paper09.py)

## 謝辞

2026年11月の原稿作成にあたり、記号的な導出、実装、日本語と英語での記述、および図を決定論的に生成するコードの作成を支援した Claude Opus 5.5、Gemini 4 Neon、GPT Astra 6.2 に感謝する。計算と図は、添付の数学的定義とデータから再現できる。内容と残る誤りの責任は著者にある。

本研究は非公式の共同プロジェクト **DPPU**（**D**onut-like topology, **P**lanck-scale compactness, **P**recession dynamics, **U**niverse）の中で進めた。
