# 问题二

## 1. 模型建立

### （1）温度分布函数

记函数：

$$
T(x,t;d_2)
$$

为第二层厚度为 $d_2$ 时的温度分布函数。

---

### （2）舒适度指标构建

采用平方根函数刻画舒适度递减趋势，可体现服装厚度增加后人体舒适性的边际损失递增特征，即随着厚度不断增加，服装笨重感与行动受限程度会加速上升。

定义归一化舒适度指数：

$$
C=C(d_2)=1-\sqrt{
\frac{
d_2-d_{2,\min}
}{
d_{2,\max}-d_{2,\min}
}
}
$$

其中：

$$
C\in[0,1]
$$

当 $d_2$ 越小时：

- $C$ 越大
- 舒适程度越高

---

### （3）隔热性能指标构建

引入热阻概念。

热阻指热量在物体中传输时，在热流路径上受到的阻碍程度，其单位为：

- K/W
- 或 ℃/W

热阻定义为：

$$
\gamma=\frac{d_n}{Ak_n}
$$

其中：

- $\gamma$：热阻
- $d_n$：第 $n$ 层厚度
- $k_n$：第 $n$ 层热导率
- $A$：横截面积

由于各层面积相同，令：

$$
A=1
$$

则热阻越大，其隔热性能越强。

定义归一化隔热性能指数：

$$
R=R(d_2)=
\frac{
\gamma_2-\gamma_{2,\min}
}{
\gamma_{2,\max}-\gamma_{2,\min}
}
$$

其中：

$$
R\in[0,1]
$$

当 $d_2$ 越大时：

- $\gamma$ 越大
- $R$ 越大
- 隔热性能越强

---

### （4）效用函数构建

采用 Cobb-Douglas 效用函数，以体现舒适性与隔热性能之间的替代关系，同时满足边际效用递减特征。

定义效用函数：

$$
U=U(d_2)=C^\alpha R^{1-\alpha}
$$

其中：

$$
U\in[0,1]
$$

---

## 2. 优化目标

求解：

$$
\max U
$$

---

## 3. 约束条件

设：

$$
L=\sum_{i=1}^{4}x_i
$$

则约束条件包括：

### （1）工作时段内皮肤外侧最高温度不超过 47℃

$$
\max_{0\le t\le60\text{min}}
T(L,t;d_2)
\le47^\circ C
$$

---

### （2）超过 44℃ 的累计时间不超过 5min

由于皮肤外侧温度关于时间单调不减，因此定义：

$$
t_{44}=\min\{t\mid T(L,t;d_2)\ge44^\circ C\}
$$

则超过 $44^\circ C$ 的累计时间为：

$$
\tau_{44}=60-t_{44}
$$

约束条件可写为：

$$
\tau_{44}\le5\text{min}
$$

等价地：

$$
t_{44}\ge55\text{min}
$$

---

### （3）厚度范围约束

$$
0.6\text{mm}\le d_2\le25\text{mm}
$$

---

## 4. 临界厚度求解

由问题一可知：

- 皮肤外侧温度关于时间单调不减
- 温度关于第二层厚度单调递减

因此可利用二分法或黄金分割法求解满足安全条件的最小临界厚度。

设：

$$
\widetilde d_{2,\min}
$$

为满足全部约束条件的最小临界厚度。

---

### （1）最高温度约束对应临界值

求解：

$$
\begin{cases}
T(L,60\text{min};D_1)=47^\circ C\\
0.6\text{mm}\le d_2\le25\text{mm}
\end{cases}
$$

则：

$$
d_2\ge D_1
$$

---

### （2）超温时间约束对应临界值

求解：

$$
\begin{cases}
T(L,55\text{min};D_2)=44^\circ C\\
0.6\text{mm}\le d_2\le25\text{mm}
\end{cases}
$$

则：

$$
d_2\ge D_2
$$

---

最终安全可行域下界为：

$$
\widetilde d_{2,\min}=\max\{D_1,D_2\}
$$

---

## 5. 综合优化模型

单纯追求最小厚度虽然能够满足安全条件，但隔热性能可能较差。

因此构建综合效用函数，在安全可行域内同时考虑：

- 舒适度
- 隔热性能

从而求得综合意义下的最优厚度。

约束区域变为：

$$
\widetilde d_{2,\min}
\le d_2\le25\text{mm}
$$

综合模型如下：

$$
\begin{aligned}
&\max U\\
&s.t.
\quad
\widetilde d_{2,\min}
\le d_2\le25\text{mm}
\end{aligned}
$$

其中：

$$
\begin{cases}
C=
1-\sqrt{
\dfrac{
d_2-d_{2,\min}
}{
d_{2,\max}-d_{2,\min}
}
}\\
\\
R=
\dfrac{
\gamma_2-\gamma_{2,\min}
}{
\gamma_{2,\max}-\gamma_{2,\min}
}
\\
\\
U=C^\alpha R^{1-\alpha}
\end{cases}
$$

最终利用：

- 二分法
- 黄金分割法

即可求解：

$$
(U_{\max},d_2^\ast)
$$

---

# 问题三

## 1. 模型建立

记：

$$
T(x,t;d_2,d_4)
$$

为第二层厚度为 $d_2$、第四层厚度为 $d_4$ 时的温度分布函数。

---

### （1）舒适度指标

定义：

$$
C=
1-\sqrt{
\frac{
(d_2+d_4)-(d_{2,\min}+d_{4,\min})
}{
(d_{2,\max}+d_{4,\max})-(d_{2,\min}+d_{4,\min})
}
}
$$

其中：

$$
C\in[0,1]
$$

总厚度越小：

- 舒适度越高

---

### （2）隔热性能指标

定义总热阻：

$$
\gamma_n=\frac{d_n}{Ak_n}
$$


且热阻具有可加性，故此时：

$$
\\gamma_{\text{sum}} = \frac{d_2}{A k_2} + \frac{d_4}{A k_4}\
$$

采用归一化形式定义综合隔热性能指标：

$$
R=
\frac{
(\gamma_2+\gamma_4)-(\gamma_{2,\min}+\gamma_{4,\min})
}{
(\gamma_{2,\max}+\gamma_{4,\max})-(\gamma_{2,\min}+\gamma_{4,\min})
}
$$

其中：

$$
R\in[0,1]
$$

---

### （3）效用函数

采用 Cobb-Douglas 效用函数：

$$
U=C^\alpha R^{1-\alpha}
$$

---

## 2. 优化目标

$$
\max U
$$

---

## 3. 约束条件

### （1）工作时段内最高温度不超过 47℃

$$
\max_{0\le t\le30\text{min}}
T(L,t;d_2,d_4)
\le47^\circ C
$$

---

### （2）超过 44℃ 的累计时间不超过 5min

定义：

$$
t_{44}=\min\{t\mid T(L,t;d_2,d_4)\ge44^\circ C\}
$$

则：

$$
\tau_{44}=30-t_{44}
$$

约束条件为：

$$
\tau_{44}\le5\text{min}
$$

等价地：

$$
t_{44}\ge25\text{min}
$$

---

### （3）厚度约束

$$
0.6\text{mm}\le d_2\le25\text{mm}
$$

$$
0.6\text{mm}\le d_4\le6.4\text{mm}
$$

---

## 4. 综合优化模型

综合模型如下：

$$
\begin{aligned}
&\max U\\
&s.t.
\begin{cases}
0.6\text{mm}\le d_2\le25\text{mm}\\
0.6\text{mm}\le d_4\le6.4\text{mm}
\end{cases}
\end{aligned}
$$

其中：

$$
\begin{cases}
C=
1-\sqrt{
\dfrac{
(d_2+d_4)-(d_{2,\min}+d_{4,\min})
}{
(d_{2,\max}+d_{4,\max})-(d_{2,\min}+d_{4,\min})
}
}
\\
\\
R=
\dfrac{
(\gamma_2+\gamma_4)-(\gamma_{2,\min}+\gamma_{4,\min})
}{
(\gamma_{2,\max}+\gamma_{4,\max})-(\gamma_{2,\min}+\gamma_{4,\min})
}
\\
\\
U=C^\alpha R^{1-\alpha}
\end{cases}
$$

---

## 5. 求解方法

采用 SLSQP 算法进行求解。

SLSQP 本质属于梯度型约束优化算法，其核心思想为：

- 将原非线性问题局部近似为二次规划子问题
- 利用拉格朗日函数与约束梯度信息不断修正搜索方向
- 最终逐步逼近最优解

由于本文目标函数与约束条件均连续可导，因此 SLSQP 能够有效利用梯度信息提高求解效率与收敛稳定性。
