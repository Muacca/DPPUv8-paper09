# 付録E　半直線上での閾値の処方

角方向の固有値 $`\lambda`$ を固定する。無限遠で $`q=r+O(r^{-3})`$ となるように固有距離座標の原点を選ぶ。動径拘束と[式(8)](DPPUv8-paper09_sec04.md#eq8)から

```math
b_s=\lambda/r+O_\lambda(r^{-3}),\qquad
 b_s^2-b_s'=\lambda(\lambda+1)r^{-2}+U_s(r),\qquad U_s=O_\lambda(r^{-4}).
```

である。エネルギー $`-\epsilon^2<0`$ では、逆二乗ポテンシャルの基準解に対して減衰の漸近係数をそろえた解を選ぶ。 $`\nu=|\lambda+1/2|\ge1`$ に対して減衰する基準解は

```math
\Phi_\epsilon(r)=\frac{\epsilon^\nu\sqrt r\,K_\nu(\epsilon r)}
                         {2^{\nu-1}\Gamma(\nu)} .
```

である。 $`\epsilon\downarrow0`$ でこの解は各点において $`r^{1/2-\nu}`$ に近づく。 $`\Phi_\epsilon`$ に対する解の比に定数変化法を適用すると、 $`y\ge r`$ で Volterra 核の絶対値は次のように抑えられる。

```math
|U_s(y)|\,\Phi_\epsilon(y)^2\int_r^y\Phi_\epsilon(v)^{-2}\mathrm d v
 \ \le\ |U_s(y)|\,y/(2\nu).
```

実際、 $`x^\nu K_\nu(x)`$ が減少することから $`\Phi_\epsilon(y)^2/\Phi_\epsilon(v)^2\le(y/v)^{1-2\nu}`$ である。可積分な優関数があるので Volterra 級数には優収束定理を適用できる。固定したコンパクトな半径区間では、積分の下端を微分し、初期値問題の解の連続依存性を用いると一階導関数も制御できる。この証明には $`\nu=1`$ となる負の最低殻も含まれるが、角運動量を無限大に送ったときの一様な評価は与えない。

$`\lambda>0`$ のとき、零エネルギーで漸近的に劣位の解は $`e^{-B_s}`$ で、 $`B_s=\lambda t-s\int_0^tX_E\mathrm d t`$ とする。その漸近係数を規格化すると、両枝に共通の因子を除いて $`e^{-sJ_\infty}`$ が掛かる。喉における Dirichlet 条件の Jost 値の枝間比は $`e^{-2J_\infty}`$ となる。 $`\lambda<0`$ のときは解 $`e^{-B_s(\tau)}\int_\tau^\infty e^{2B_s(v)}\mathrm d v`$ を用いる。その漸近係数は $`e^{-sJ_\infty}`$ に比例するので、規格化による因子は $`e^{sJ_\infty}`$ である。喉での Robin 微分 $`A_s f(0)`$ は両枝で共通の非零の符号をもち、その比は $`e^{2J_\infty}`$ である。ここで対数比の角方向 Abel 和を取ると $`-2J_\infty\eta_0`$ を得る。この座標に依存しない比較では、殻ごとに外部での規格化を二つの枝で共通に固定した。全殻にわたる極限と Jost 閾値の極限の順序を入れ替える操作は用いていない。

[英語原文](LaTeX/appendices/appE.tex)

[← 付録D](DPPUv8-paper09_appD.md)　[目次](README_ja.md)　[付録F →](DPPUv8-paper09_appF.md)
