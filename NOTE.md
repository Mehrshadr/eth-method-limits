# How far can shared linear summaries go? Lower bounds for the method of "The Exponential Time Hypothesis is False"

**Authors:** Mehrshad (with Claude)
**Date:** 30 October 2026

**Status.** This note has not been peer reviewed. Its arguments were developed with an AI assistant (Claude). Earlier versions of Theorem A and Theorem C were each checked by a separate AI reviewer, but this written note as a whole has not been independently reviewed. Lemma 2.1 (the triangular fooling-set lemma) is machine-checked in Lean 4 with Mathlib. That is weaker than expert human review, so read it as a careful preprint draft, not as an established result. Every statement carries one of these labels:

- **[R] rigorous:** a complete proof is given here.
- **[R, c]:** rigorous, but the constants have not been optimised.
- **[S] sketch:** the argument is outlined and some steps are not written out.
- **[E] experimental:** supported only by computer searches. There is no proof.

---

## Abstract

OpenAI's paper "The Exponential Time Hypothesis is False" (17 Oct 2026, revised 25 Oct 2026; called *the paper* below) counts the satisfying assignments of a bounded-occurrence 3-CNF in the following way:

- It replaces each clause input by a low-degree polynomial family of linear forms $P_o(T)$.
- At $Q$ evaluation points it makes the variable groups interact only through a few shared integer linear measurements, and it enumerates every value those measurements can take.
- It recovers the count at $T=0$ from Taylor data by inclusion–exclusion.

We fix the paper's architecture (Sections 2.4–2.7, Definition 3.1, Lemma 3.2, Lemma 3.3) and prove lower bounds that hold for **every** choice of its free parameters. The parameters are the schedules, the number of groups $r$, $Q$, $R$, the "leak" matrix $M_0$ and the measurements. On bounded-occurrence 3-CNFs whose variable graph is an edge expander, we show:

- **(A)** The one-level routine needs $2^{\Omega(\sqrt n)}$ time even with measurements that vary from point to point.
- **(B, C)** The number of measurements is at least $n\cdot |X|^{-O(R/Q)}$, where $X$ is the set of schedule points.
- **(D)** The number of measurements is at least $n\cdot e^{-O(R^2/Q)}$. So the paper's "bell curve" relation $R\gtrsim\sqrt{Q\log(1/\tau)}$ is forced.
- **(E)** Combining these with the Taylor cost, every design in the framework costs $2^{\Omega(n(\log\log n)^2/\log n)}$. With $Q=O(\log n)$, which covers the paper's own parameter range, it costs $2^{\Omega(n/\sqrt{\log n})}$.

None of this contradicts the paper. Its Section 2–3 bound is $2^{O(n(\log\log n)^2/\sqrt{\log n})}$, and Appendix D uses a different, recursive method that lies outside the framework. Section 8 looks at four ways out of the framework. Section 9 summarises computer experiments; they are evidence, not proofs.

---

## 1. The framework

### 1.1 Instances

Let $F$ be a CNF on the variables $[n]$ with clause set $\mathcal F$ and $m=|\mathcal F|$. Every clause has at most 3 literals, and every variable appears in at most $B$ clauses. Section 2.1 of the paper has $B=3$; we keep $B$ as a symbol. Let $I_e$ be the variable set of clause $e$. The **occurrences** are
$$\mathcal O=\{(e,v):v\in I_e\},\qquad |\mathcal O|\le Bn.$$
The **variable graph** $G_F$ has vertex set $[n]$, and $v\sim w$ when some clause contains both. Each vertex has degree at most $2B$. $G_F$ has **edge expansion** $\varphi>0$ if
$$e(U,[n]\setminus U)\ \ge\ \varphi\,|U|\qquad\text{for every } U\subseteq[n] \text{ with } |U|\le n/2 .$$
Such families exist with $B=3$. Take a 3-regular graph with edge expansion bounded away from 0 (random cubic graphs have this property with high probability) and add the clause $(x_u\vee x_w)$ for each edge $uw$.

### 1.2 The architecture (paper, Sections 2.3–2.7 and 3.1–3.5)

**(F1) Points and schedules.** Fix integers $1\le R<Q$ and $r\ge 2$. Each variable $v$ gets a home schedule $a_v\in[r]^Q$, and each clause $e$ gets a reader schedule $d_e\in[r]^Q$. At point $t$ the groups are $S_{t,g}=\{v:a_v(t)=g\}$, and clause $e$ is assigned to group $\sigma_e(t)=d_e(t)$ (Definition 2.1(1)–(2), Section 2.7 Step 4). The routine scans the cube of each group, so $s_t:=\max_g|S_{t,g}|$ costs a factor $2^{s_t}$. The paper requires $|S_{t,g}|\le 2n/r$. Write $\operatorname{dist}$ for Hamming distance on $[r]^Q$. The **blocked set** of an occurrence $o=(e,v)$ is
$$E_o=\{t:a_v(t)\ne d_e(t)\},\qquad |E_o|=\operatorname{dist}(a_v,d_e).$$

**(F2) Coordinates.** There is a finite coordinate set $\mathcal X$. Each coordinate $u$ has its own home and reader schedules $a_u,d_u\in[r]^Q$ and a fill $\mathrm{var}(u)\in[n]\cup\{0\}$. Each occurrence $o=(e,v)$ owns a coordinate $u(o)$ with $\mathrm{var}=v$, $a_{u(o)}=a_v$ and $d_{u(o)}=d_e$. The filled vector is $\zeta(x)=Fx$, where $F[u,v]=1$ exactly when $\mathrm{var}(u)=v$. In the paper, $\mathcal X$ is $h$ copies of $B=([r]\times[r])^Q$ and unused coordinates are filled with 0 (Section 2.7, Steps 1 and 4).

**(F3) Leak matrix.** $M_0\in\mathbb Q^{\mathcal X\times\mathcal X}$ satisfies the support condition
$$\textbf{(S)}\qquad M_0[u,w]\ne 0\ \Longrightarrow\ \operatorname{dist}(a_w,d_u)\le R-1 .$$
This is the second condition of (2.14). The replacement row of coordinate $u$ is (2.16):
$$\langle P_u(T),\zeta\rangle=((I-M_0)\zeta)_u+\sum_{w}M_0[u,w]\,\phi_{u,w}(T)\,\zeta_w .$$
Here $\phi_{u,w}(0)=1$, $\deg\phi_{u,w}\le R-1$, and $\phi_{u,w}$ vanishes on $E_{u,w}=\{t:a_w(t)\neq d_u(t)\}$. Such a $\phi$ exists exactly when $|E_{u,w}|\le R-1$, so (S) is forced and not a choice. The occurrence row is $P_o(T)=F^{\mathsf T}P_{u(o)}(T)$, which satisfies $P_o(0)=e_v$ and $\deg P_o\le R-1\le R$.

**(F4) Measurements.** There is one integer matrix $C$ on $\mathcal X$ whose row space contains every row of $I-M_0$ (first condition of (2.14)). At each $t$, the measurements of $C\zeta(x)$ are split into group pieces $C_{t,g}$ (Section 2.7 Step 4, Definition 2.1(3)). So the number of measurement rows is $\ge\operatorname{rank}C\ge\operatorname{rank}(I-M_0)$.

**(F5) Separation.** At each $t$ and for each inclusion–exclusion pair, the routine of Section 2.4 / Lemma 3.2 enumerates measurement tuples $\vec y$. For each tuple it scans every group cube.

**(F6) Recovery.** $q=\lfloor 3nR/Q\rfloor$ and $J=q+1$ (Section 2.5, (3.5)). At each point, the routine computes Taylor data through Lemma 3.3, which uses all pairs $U\subseteq E\subseteq[n]$ with $|E|\le q$. There are $C_q=\sum_{s\le\min(q,n)}2^s\binom ns$ such pairs, as in (3.14).

The paper's construction is one instance: $M_0=V_0f(W_0)$ with (2.15) giving (S), and $C$ made of the exceptional eigenvectors.

**Two variants are deliberately excluded from (F1)–(F6).** (i) Definition 3.1 allows the measured subspace $C_{t,g}$ to change with $t$ in an arbitrary way and does not ask for one global $M_0$. Theorem A covers that generality; Theorems B–D and the corollaries do **not**. (ii) Any recovery other than (F6) is outside the framework (Section 8).

### 1.3 Cost accounting

**Lemma 1.1 [R].** Let $L:\mathbb Q^N\to\mathbb Q^k$ be linear of rank $\rho$. Then $|L(\{0,1\}^N)|\ge 2^\rho$.

*Proof.* Choose $\rho$ linearly independent columns $L e_{i_1},\dots,L e_{i_\rho}$. On the subcube where only $x_{i_1},\dots,x_{i_\rho}$ vary and every other coordinate is 0, $L$ is injective. So it takes $2^\rho$ distinct values there. $\square$

**Lemma 1.2 [R] (variable-level leak).** Define $N\in\mathbb Q^{\mathcal O\times[n]}$ by $N[o,w]=\sum_{c:\mathrm{var}(c)=w}M_0[u(o),c]$. Let $E[o,w]=1$ exactly when $w=v(o)$. Then:

1. $N[o,w]\ne0$ implies $\operatorname{dist}(a_w,d_{e(o)})\le R-1$.
2. $\operatorname{rank}(CF)\ge\operatorname{rank}(E-N)$.

*Proof.* (1) A nonzero term needs some $c$ with $\mathrm{var}(c)=w$ and $M_0[u(o),c]\neq0$. Every coordinate holding $w$ has home schedule $a_w$, so (S) gives the claim. (2) Write $I-M_0=KC$. Then $(I-M_0)F=KCF$, so $\operatorname{rank}(CF)\ge\operatorname{rank}((I-M_0)F)$. Row $u(o)$ of $F$ is $e_{v(o)}^{\mathsf T}$, so row $u(o)$ of $(I-M_0)F$ is $(E-N)[o,\cdot]$. $\square$

**Work.** Let $s=\max_t s_t$. In (F5), the literal enumeration of Lemma 3.2 visits $(2D+1)^{\#\text{rows}}\ge 3^{\operatorname{rank}(I-M_0)}$ tuples. An enumerator that visits only *realised* tuples still visits at least $|\{C\zeta(x)\}|\ge2^{\operatorname{rank}(CF)}\ge 2^{\operatorname{rank}(E-N)}$ of them, by Lemmas 1.1 and 1.2. For each realised tuple it scans at least $2^{s_t}$ vectors. Hence the work is at least
$$\mathrm{Work}\ \ge\ C_q\cdot 2^{\operatorname{rank}(E-N)}\cdot 2^{s}. \tag{1.1}$$
If some group has more than $n/2$ variables, the work is $\ge 2^{n/2}$. **So from now on we assume every group has at most $n/2$ variables.**

Two facts about $C_q$ [R]:
- $C_q\ge 2^{\min(q,n)}$: keep only the term $s=\min(q,n)$.
- $C_q\ge\binom n{q'}\ge (n/q')^{q'}$ for every $q'\le\min(q,n)$.

---

## 2. The triangular fooling-set lemma

Let $M_0\in\mathbb F^{\mathcal X\times\mathcal X}$ over any field $\mathbb F$, with support inside a digraph $\Gamma$ on $\mathcal X$: $M_0[u,w]\ne0$ implies $u\to w$ in $\Gamma$. A loop is an arc $u\to u$.

**Lemma 2.1 [R].** Let $S\subseteq\mathcal X$. Suppose no vertex of $S$ has a loop and the induced digraph $\Gamma[S]$ is acyclic. Then $\operatorname{rank}(I-M_0)\ge|S|$.

*Proof.* Order $S$ topologically, so every arc of $\Gamma[S]$ goes from an earlier vertex to a later one. In that order, the principal submatrix $M_0[S,S]$ has zero diagonal (no loops) and no entries below the diagonal (arcs go forward). It is strictly upper triangular. So $(I-M_0)[S,S]$ is upper unitriangular, with determinant 1 and rank $|S|$. A submatrix has rank at most the rank of the whole matrix. $\square$

The same lemma applies to $E-N$ (Lemma 1.2) once we fix a set $S$ of variables and one occurrence $o(v)$ for each $v\in S$. The submatrix $(E-N)[o(S),S]$ equals $I-N'$ with $N'[v,w]=N[o(v),w]$. It has an arc $v\to w$ only if $\operatorname{dist}(a_w,d_{e(o(v))})\le R-1$. We write $\Gamma$ for this relation in both settings: $o\to o'$ when $\operatorname{dist}(a_{v(o')},d_{e(o)})\le R-1$.

---

## 3. Theorem A: the one-level barrier

**Theorem A [R, c].** Let $F$ be a 3-CNF on $n$ variables in which each variable appears in at most 3 clauses, and whose variable graph has edge expansion $\varphi>0$. Take any data satisfying Definition 3.1, with **arbitrary** $t$-dependent $C_{t,g}$ and $q=\lfloor 3nR/Q\rfloor$. Then the routine of Proposition 3.4 examines at least
$$2^{\min\{\sqrt{\varphi n/6},\ \lfloor \varphi n/8\rfloor,\ n/2\}}=2^{\Omega(\sqrt n)}$$
Boolean vectors. This holds even if each call of Lemma 3.2 enumerates only realised measurement tuples. With literal enumeration, $3$ replaces $2$ as the base of the measurement factor below.

*Proof.* Put $\varepsilon=\varphi/24$. If some group ever has more than $n/2$ variables, scanning its cube already costs $2^{n/2}$, so assume every group has at most $n/2$ variables.

*Case $R\ge\varepsilon Q$.* Then $q\ge\lfloor 3\varepsilon n\rfloor=\lfloor\varphi n/8\rfloor$, and $C_q\ge 2^{\min(q,n)}$.

*Case $R<\varepsilon Q$.* For an occurrence $o=(e,v)$, let $p_o(T)$ be the coefficient of $X_v$ in $P_{e,v}(T)$. It has degree $\le R$ and $p_o(0)=1$, so it is a nonzero polynomial and its zero set $Z_o\subseteq[Q]$ has $|Z_o|\le R$. Then
$$\sum_t\#\{o:t\in Z_o\}=\sum_o|Z_o|\le R|\mathcal O|<\varepsilon Q|\mathcal O|.$$
So some point $t$ has fewer than $\varepsilon|\mathcal O|\le 3\varepsilon n=\varphi n/8$ occurrences with $p_o(t)=0$. Fix that $t$.

1. *Cut edges.* Every group has at most $n/2$ variables, so
   $$\sum_g e(S_{t,g},\overline{S_{t,g}})\ge\varphi\sum_g|S_{t,g}|=\varphi n .$$
   Each cut edge is counted twice in this sum, so there are at least $\varphi n/2$ cut edges.

2. *Blocked occurrences.* A cut edge $\{v,w\}$ lies in some clause $e$. Since $a_v(t)\ne a_w(t)$, at least one of $(e,v)$ and $(e,w)$ is *blocked* at $t$, meaning its variable is not in group $\sigma_e(t)$. Map the edge to such an occurrence. An occurrence $(e,v)$ receives at most two edges, one for each other variable of $e$. So at least $\varphi n/4$ occurrences are blocked at $t$, and more than $\varphi n/4-\varphi n/8=\varphi n/8$ of them also have $p_o(t)\ne0$.

3. *Touched groups.* Take such an occurrence, with $h=a_v(t)\ne\sigma_e(t)$. Definition 3.1(4) puts $P_{e,v}(t)|_{S_{t,h}}$ in the row space of $C_{t,h}$. That restricted row has the nonzero entry $p_o(t)$ at $v$. So $C_{t,h}\ne0$ and group $h$ is *touched*. A group with at most $s$ variables holds at most $3s$ occurrences. So more than $\varphi n/(24s)$ groups are touched.

4. *Cost.* On each touched group, $C_{t,h}$ is a nonzero linear map. By Lemma 1.1 it takes at least 2 values on $\{0,1\}^{S_{t,h}}$. The realised tuples form a product over groups, so there are at least $2^{\varphi n/(24s)}$ of them. For each realised tuple, Lemma 3.2 scans all of $\{0,1\}^{S_{t,g}}$ for every $g$, which is at least $2^s$ vectors. The work is therefore at least $2^{s+\varphi n/(24s)}\ge 2^{2\sqrt{\varphi n/24}}=2^{\sqrt{\varphi n/6}}$ by AM–GM. $\square$

**Remarks.**
1. The brief stated "touched groups $\ge(\alpha-3\varepsilon)m/(3s)$ with $\varepsilon<\alpha/3$". With the normalisation $e(U,\bar U)\ge\varphi|U|$ used here, the correct count is $(\varphi/4-3\varepsilon)n/(3s)$, which needs $\varepsilon<\varphi/12$. We use $\varepsilon=\varphi/24$.
2. Theorem A is far below the paper's upper bounds. Its value is generality: it allows arbitrary $t$-dependent measurements and still rules out, for example, $2^{n^{0.49}}$ for the one-level routine.
3. **Appendix D escapes it.** Its recursion (Theorem D.1, Section D.4) contracts each piece recursively instead of scanning a group's cube by brute force. The $2^{s}$ factor, and with it the $\sqrt n$ balance, disappears.

---

## 4. Theorem B: the two-room ceiling

This is a warm-up that shows the mechanism behind Theorems C and D. Take $r=2$ and identify the rooms $\{A(t),B(t)\}$ with $\{0,1\}$, so schedules are vectors in $\{0,1\}^Q$. Let $U\subseteq\mathcal X$ be a set of $N$ coordinates that are **fully blocked**: $d_u=\mathbf 1\oplus a_u$, so the reader is in the other room at every point. By (S), an arc $u\to w$ inside $U$ requires
$$\operatorname{dist}(a_w,\mathbf 1\oplus a_u)\le R-1 ,$$
that is, $a_w$ is a *near-antipode* of $a_u$. This relation is symmetric. Let $H$ be the undirected graph on $U$ with these edges. It has no loops, because $\operatorname{dist}(a_u,\mathbf 1\oplus a_u)=Q>R-1$.

**Lemma 4.1 [R].** Every odd cycle of $H$ has length at least $Q/(R-1)$. (We assume $R\ge2$. If $R=1$, $H$ has no edges.)

*Proof.* Each step of a walk replaces $a$ by $\mathbf 1\oplus a\oplus\eta_i$ with $|\eta_i|\le R-1$. After $\ell$ steps the vector is $a\oplus(\ell\bmod 2)\mathbf 1\oplus\bigoplus_i\eta_i$. A closed walk with $\ell$ odd therefore needs $\bigoplus_i\eta_i=\mathbf 1$. This forces $\sum_i|\eta_i|\ge Q$, hence $\ell(R-1)\ge Q$. $\square$

**Lemma 4.2 [R] (independence number under odd girth).** Let $G$ be a graph on $N$ vertices with no odd cycle of length at most $2k-1$, where $k\ge1$. Then $\alpha(G)\ge N^{1-1/k}/2$.

*Proof.* We use BFS ball growing. Let $G'$ be the current induced subgraph with $n'\le N$ vertices, and put $\Delta=N^{1/k}$. Pick a vertex $v$ and let $b_i=|B_i(v)|$ be the size of the ball of radius $i$ in $G'$. Since $b_0=1$ and $b_k\le n'\le N=\Delta^k$, some $i\in\{0,\dots,k-1\}$ has $b_{i+1}\le\Delta b_i$. (Otherwise $b_k>\Delta^k\ge n'$.)

For $j\le i\le k-1$, the layer $L_j$ (distance exactly $j$) is independent. An edge inside $L_j$ together with two shortest paths would form a closed walk of odd length $2j+1\le 2k-1$, and every odd closed walk contains an odd cycle of at most its length. BFS edges join only equal or consecutive layers. So the even layers of $B_i(v)$ form an independent set, and so do the odd layers. The larger of the two, $I_v$, has $|I_v|\ge b_i/2$. All its neighbours lie in $B_{i+1}(v)$.

Add $I_v$ to the output, delete $B_{i+1}(v)$, and repeat. Odd-cycle-freeness passes to induced subgraphs. Sets chosen later lie outside the deleted balls, so they are not adjacent to earlier ones, and the union is independent. Each round deletes at most $\Delta b_i\le 2\Delta|I_v|$ vertices, and the rounds delete all $N$ vertices. Hence $\sum|I_v|\ge N/(2\Delta)=N^{1-1/k}/2$. $\square$

**Theorem B [R, c].** Assume $R\ge2$ and put $g=\lceil Q/(R-1)\rceil$ and $k=\lfloor (g-1)/2\rfloor\ge1$. For every $M_0$ satisfying (S),
$$\operatorname{rank}(I-M_0)\ \ge\ \alpha(H)\ \ge\ N^{1-1/k}/2 .$$
So if at most $\tau N$ measurements cover these coordinates (that is, $\operatorname{rank}(I-M_0)\le\tau N$), then
$$N\ \ge\ (2\tau)^{-k}.$$
When $Q\ge 2R(R-1)$ we have $k\ge Q/(2(R-1))-1\ge Q/(2R)$, and this gives $N\ge(1/(2\tau))^{Q/(2R)}$ for $\tau\le1/2$.

*Proof.* An independent set of $H$ has no arcs in either direction and no loops. Apply Lemma 2.1, then Lemmas 4.1 and 4.2 with $2k+1\le g$. $\square$

**Corrections to the brief.** The form $N\ge(c/\tau)^{Q/(2R)}$ holds with $c=1/2$ only in the range $Q\ge 2R(R-1)$. In general use $(2\tau)^{-\lfloor(\lceil Q/(R-1)\rceil-1)/2\rfloor}$. If the coordinates are only $K$-nearly blocked, $\operatorname{dist}(a_u,d_u)\ge Q-K$, the same proof works with $R-1$ replaced by $R-1+K$. **Reading:** a measured fraction $\tau$ forces $\log N\gtrsim (Q/R)\log(1/\tau)$. Since $N\le|\mathcal X|$, we get $\tau\ge |\mathcal X|^{-O(R/Q)}$.

---

## 5. Theorem C: the multi-room ceiling

Let $\mu:=\varphi/(8B)$. Call an occurrence $o$ **heavy** if $|E_o|\ge\mu Q$. Let $X$ be the multiset of schedule points $\{a_v\}_{v\in[n]}\cup\{d_e\}_{e\in\mathcal F}$, with the Hamming pseudometric. Then $|X|=n+m$, which is at most $(1+B)n$.

**Lemma 5.1 [R] (many heavy occurrences).** Under (F1), with expansion $\varphi$ and all groups of size at most $n/2$, at least $\varphi n/8$ occurrences are heavy.

*Proof.* Step 2 of Theorem A's proof applies at every $t$, without the zero-set restriction: at least $\varphi n/4$ occurrences are blocked at each $t$. So $\sum_o|E_o|\ge\varphi nQ/4$. The light occurrences contribute less than $\mu Q\cdot|\mathcal O|\le\mu QBn=\varphi nQ/8$. The heavy ones therefore contribute at least $\varphi nQ/8$, each at most $Q$, so there are at least $\varphi n/8$ of them. $\square$

**Exponential-shift partition (Miller–Peng–Xu).** Fix $\beta>0$. Give every index $x\in X$ an independent shift $\delta_x\sim\mathrm{Exp}(\beta)$. Assign every point $y$ to the centre
$$\mathrm{ctr}(y)=\arg\max_{x\in X}\big(\delta_x-\operatorname{dist}(x,y)\big).$$
Ties have probability 0. Independently, pick a uniformly random order $\pi$ of $X$.

**Lemma 5.2 [R] (gap lemma).** Let $c_x$ be real constants and $Z_x=\delta_x+c_x$. Let $\mathrm{gap}$ be the largest $Z_x$ minus the second largest. Then $\Pr[\mathrm{gap}>g]\ge e^{-\beta g}$ for every $g\ge0$.

*Proof.* Fix $x$ and condition on $(\delta_{x'})_{x'\ne x}$. Put $a=\max_{x'\ne x}Z_{x'}-c_x$. Let $A_x=\{\delta_x>a\}$ (x is the maximiser) and $G_x=\{\delta_x>a+g\}$ (x wins with gap $>g$).
- If $a\ge0$, memorylessness gives $\Pr[G_x]=e^{-\beta g}\Pr[A_x]$.
- If $a<0$, then $\Pr[A_x]=1$ and $\Pr[G_x]=e^{-\beta\max(a+g,0)}\ge e^{-\beta g}$.

The events $G_x$ are disjoint and their union is $\{\mathrm{gap}>g\}$. Hence $\Pr[\mathrm{gap}>g]=\sum_x\Pr[G_x]\ge e^{-\beta g}\sum_x\Pr[A_x]=e^{-\beta g}$. $\square$

**Padding [R].** Suppose the gap at $y$ (with $c_x=-\operatorname{dist}(x,y)$) exceeds $2\rho$. Then every point $z$ with $\operatorname{dist}(y,z)\le\rho$ has $\mathrm{ctr}(z)=\mathrm{ctr}(y)$. Indeed, each score changes by at most $\rho$ when moving from $y$ to $z$. By Lemma 5.2, $\Pr[B(y,\rho)\text{ padded}]\ge e^{-2\beta\rho}$.

**Separation [R].** Suppose $a_v$ and $d_e$ have the same centre $c$. Every point $y$ is a candidate centre for itself, with score $\delta_y\ge0$. So $\delta_c-\operatorname{dist}(c,a_v)\ge0$ and $\delta_c-\operatorname{dist}(c,d_e)\ge0$, which gives $\operatorname{dist}(a_v,d_e)\le2\delta_c$. For a heavy $o=(e,v)$ this requires $\max_x\delta_x\ge\mu Q/2$, an event of probability at most $|X|e^{-\beta\mu Q/2}$.

**Theorem C [R, c].** Assume all groups have at most $n/2$ variables, the variable graph has expansion $\varphi$, and $8(R-1)\le\mu Q$. Then for every $M_0$ satisfying (S),
$$\operatorname{rank}(I-M_0)\ \ge\ \tfrac{3\varphi n}{64}\,(4|X|)^{-8(R-1)/(\mu Q)}\ \ge\ \tfrac{\varphi n}{32}\,(4|X|)^{-8R/(\mu Q)},$$
$$\operatorname{rank}(E-N)\ \ge\ \tfrac{3\varphi n}{64B}\,(4|X|)^{-8(R-1)/(\mu Q)} .$$

*Proof.* Set $\beta=4\ln(4|X|)/(\mu Q)$ and $\rho=R-1$. Call a heavy occurrence $o=(e,v)$ **good** if three things hold:
1. the ball $B(d_e,\rho)$ is padded;
2. $\mathrm{ctr}(a_v)\ne\mathrm{ctr}(d_e)$;
3. $\pi(\mathrm{ctr}(a_v))<\pi(\mathrm{ctr}(d_e))$.

*Probability of good.* We have $e^{-2\beta\rho}=(4|X|)^{-8\rho/(\mu Q)}\ge(4|X|)^{-1}$, using $8\rho\le\mu Q$. The separation failure probability is at most $|X|(4|X|)^{-2}=1/(16|X|)\le\tfrac14e^{-2\beta\rho}$. So events 1 and 2 hold together with probability at least $\tfrac34e^{-2\beta\rho}$. Given them, event 3 has probability $1/2$, because $\pi$ is independent. Altogether $\Pr[\text{good}]\ge\tfrac38e^{-2\beta\rho}$.

*Acyclicity.* Let $S$ be the set of good occurrences and order it by $\kappa(o)=\pi(\mathrm{ctr}(a_{v(o)}))$. Take an arc $o\to o'$ inside $S$. It means $\operatorname{dist}(a_{v(o')},d_{e(o)})\le\rho$. Padding then gives $\mathrm{ctr}(a_{v(o')})=\mathrm{ctr}(d_{e(o)})$. So $\kappa(o')=\pi(\mathrm{ctr}(d_{e(o)}))>\kappa(o)$. Arcs strictly increase $\kappa$, so $\Gamma[S]$ is acyclic and has no loops. Lemma 2.1 gives $\operatorname{rank}(I-M_0)\ge|S|$ for every outcome. Hence
$$\operatorname{rank}(I-M_0)\ge\mathbb E|S|\ge\tfrac{\varphi n}8\cdot\tfrac38e^{-2\beta\rho}$$
by Lemma 5.1.

*Variable level.* Keep one good occurrence per variable. That keeps at least $|S|/B$ variables, and the same ordering argument applies to $E-N$ (Lemma 1.2(1)). $\square$

**Notes.**
1. The hypothesis $8(R-1)\le\mu Q$ is **necessary** for any nontrivial bound. If $R-1\ge\max_o|E_o|$, then $N=E$ satisfies (S) and the rank can be 0. The brief stated the bound without this hypothesis.
2. The brief's constant $\varphi n/32$ is correct for the coordinate-level rank, which counts the literal measurement rows. For *realised* measurement values one must divide by $B$.
3. $|X|\le n+\#\text{clauses}$.

---

## 6. Theorem D: the bell curve is unavoidable

Replace the clusters by a random potential. Let the $\chi_t(a)$, for $t\in[Q]$ and $a\in[r]$, be independent uniform signs $\pm1$, and define
$$\Psi(x)=\sum_t\chi_t(x_t)\qquad\text{for } x\in[r]^Q .$$
(The brief calls this $\varphi$; we write $\Psi$ because $\varphi$ is already the expansion.) If $\operatorname{dist}(x,y)\le R-1$, then $|\Psi(x)-\Psi(y)|\le 2(R-1)$.

**Lemma 6.1 [R, c] (reverse tail).** Let $N$ be even, $h\ge0$, and $h+\sqrt N\le N/4$. Then
$$\Pr[\mathrm{Bin}(N,\tfrac12)\ge N/2+h]\ \ge\ \tfrac1{\sqrt2}\,e^{-16}\,e^{-16h^2/N}.$$

*Proof.* There are three steps.
1. *Central term.* $\binom{N}{N/2}2^{-N}\ge(2N)^{-1/2}$. Write $N=2k$ and $a_k=\prod_{i\le k}\frac{2i-1}{2i}$. We show $a_k\ge 1/(2\sqrt k)$ by induction. The step needs $(2k+1)^2(k+1)\ge(2k+2)^2k$, which reduces to $4k^2+4k+1\ge4k^2+4k$.
2. *Ratio.* For $0\le m\le N/4$, put $x=2m/N\le1/2$. Then
   $$\frac{\binom N{N/2+m}}{\binom N{N/2}}=\prod_{i=1}^m\frac{N/2-i+1}{N/2+i}\ge\Big(\frac{1-x}{1+x}\Big)^m .$$
   Also $\ln\frac{1-x}{1+x}\ge-\frac{x}{1-x}-x\ge-4x$. So the ratio is at least $e^{-8m^2/N}$.
3. *Sum.* Sum over $m=h,\dots,h+L$ with $L=\lfloor\sqrt N\rfloor$. There are $L+1\ge\sqrt N$ terms, each at least $(2N)^{-1/2}e^{-8(h+L)^2/N}$. Finally $(h+L)^2\le2h^2+2N$. $\square$

**Theorem D [R, c].** Assume all groups have at most $n/2$ variables, the variable graph has expansion $\varphi$, $R\le\mu Q/4$ and $\mu Q\ge32$. Then for every $M_0$ satisfying (S),
$$\operatorname{rank}(I-M_0)\ \ge\ \tfrac{\varphi n}{8}\cdot c\cdot e^{-8R^2/(\mu Q)},\qquad c=e^{-16}/\sqrt2 .$$
The same holds for $\operatorname{rank}(E-N)$ after dividing by $B$.

*Proof.* Let $S$ be the set of heavy $o=(e,v)$ with $\Psi(d_e)-\Psi(a_v)>2(R-1)$.

*Acyclicity.* Take an arc $o\to o'=(e',w)$ inside $S$. It means $\operatorname{dist}(a_w,d_e)\le R-1$. Then
$$\Psi(a_w)\ \ge\ \Psi(d_e)-2(R-1)\ >\ \Psi(a_v).$$
So arcs strictly increase $\Psi(a_{v(\cdot)})$, and $\Gamma[S]$ has no cycles and no loops. Lemma 2.1 gives $\operatorname{rank}(I-M_0)\ge|S|$.

*Probability.* The terms with $t\notin E_o$ cancel. For $t\in E_o$, the two signs involved are distinct and independent, and different $t$ are independent. Each difference $\chi_t(d_e(t))-\chi_t(a_v(t))$ equals $2(\eta+\eta'-1)$ with $\eta,\eta'$ independent Bernoulli(1/2). So with $b=|E_o|$,
$$\Psi(d_e)-\Psi(a_v)\ \overset{d}{=}\ 2\big(\mathrm{Bin}(2b,\tfrac12)-b\big),$$
and the event is $\mathrm{Bin}(2b,\frac12)-b\ge R$.

*Applying Lemma 6.1.* Use $N=2b\ge2\mu Q$ and $h=R$. The map $b\mapsto b/2-\sqrt{2b}$ is increasing for $b\ge2$. The hypotheses give $R+\sqrt{2\mu Q}\le\mu Q/4+\mu Q/4$, using $\mu Q\ge32$. So the condition $h+\sqrt N\le N/4$ holds for every heavy $o$, and $\Pr[o\in S]\ge c\,e^{-16R^2/(2b)}\ge c\,e^{-8R^2/(\mu Q)}$. Take expectations and use Lemma 5.1. $\square$

**Consequence.** Suppose there are at most $2\tau n$ measurement rows, as in Definition 2.1(3). Then
$$R^2\ \ge\ \tfrac{\mu}{8}\,Q\ln\frac{c\varphi}{16\tau}.$$
That is, $R\gtrsim\sqrt{Q\log(1/\tau)}$. The paper's choice in Section 2.7 Step 5 ($2e^{-R^2/(8Q)}\le\tau/(3r)$) is therefore optimal up to constants: no $M_0$ under (S) beats the binomial window by more than a constant factor in the exponent.

Theorem D has no $|X|$; Theorem C has no $R^2$. Together they give $\operatorname{rank}\gtrsim n\exp\!\big(-O(\min\{R\log|X|,R^2\})/(\mu Q)\big)$.

---

## 7. Corollaries: the framework's floor

**Corollary E [R, c].** Fix $B$ and $\varphi>0$. Every design in (F1)–(F6) on an $n$-variable 3-CNF with at most $B$ occurrences per variable and variable-graph expansion $\varphi$ does work at least $2^{c\,n(\log\log n)^2/\log n}$ for all $n\ge n_0$. Here $c$ and $n_0$ depend only on $\varphi$ and $B$. In particular the exponent is $\Omega(n\log\log n/\log n)$.

*Proof.* By (1.1), $\log_2\mathrm{Work}\ge\max(\log_2C_q,\ \operatorname{rank}(E-N))$, and we may assume all groups have at most $n/2$ variables. Put $L=\ln(4|X|)\le 2\ln n$ for large $n$.

1. *Case $8(R-1)>\mu Q$.* Then $R/Q>\mu/8$, so $q\ge\lfloor 3\mu n/8\rfloor$ and $\log_2 C_q\ge\min(q,n)=\Omega(n)$.
2. *Case $8(R-1)\le\mu Q$ and $8(R-1)L/(\mu Q)\le\frac12\ln\ln n$.* Theorem C gives $\operatorname{rank}(E-N)\ge\frac{3\varphi}{64B}\,n/\sqrt{\ln n}$.
3. *Otherwise.* $R/Q>\mu\ln\ln n/(32\ln n)$, so $q\ge q_0-1$ with $q_0=\frac{3\mu}{32}\,n\ln\ln n/\ln n$. Take $q'=\min(q,\lfloor n/2\rfloor)$. The function $x\ln(n/x)$ increases on $(0,n/e]$, so its minimum on $[q_0-1,n/2]$ is attained at an endpoint. Hence
   $$\ln C_q\ge q'\ln(n/q')\ge\min\{(q_0-1)\ln\tfrac n{q_0-1},\ \tfrac n2\ln2\}=\Omega\big(n(\ln\ln n)^2/\ln n\big).$$

The smallest of the three cases is case 3. $\square$

**Corollary F [R, c] (short schedules).** Suppose in addition that $Q\le K\log n$. This holds in the paper, where $D=r^{2Q}\le(2n)^{1/128}$ (Section 2.7 Step 5). Then the work is at least $2^{c'n/\sqrt{\log n}}$.

*Proof.* Use Theorem D in place of Theorem C.
- If $R>\mu Q/4$, or if $\mu Q<32$ (so $R/Q\ge1/Q>\mu/32$), then $q=\Omega(n)$.
- Otherwise, either $8R^2/(\mu Q)\le\frac12\ln\ln n$, and the rank is $\Omega(n/\sqrt{\ln n})$;
- or $R/Q>\sqrt{\mu\ln\ln n/(16K\ln n)}$, and $\log_2C_q\ge\min(q,n)=\Omega(n\sqrt{\ln\ln n/\ln n})$. $\square$

**Contrast with the paper.**

| Bound | Exponent |
|---|---|
| Paper, Sections 2–3 (Section 2.6, (2.12)), achieved | $O(n(\log\log n)^2/\sqrt{\log n})$ |
| Paper, Appendix D, fixed density, achieved (outside the framework) | $O(n/\sqrt{\log n})$ |
| This note, Corollary E, any design in the framework | $\Omega(n(\log\log n)^2/\log n)$ |
| This note, Corollary F, any design in the framework with $Q=O(\log n)$ | $\Omega(n/\sqrt{\log n})$ |

So within its own range of $Q$, the paper's one-level design is optimal up to a $(\log\log n)^2$ factor. Going below $n/\sqrt{\log n}$ inside the framework requires $Q\gg\log n$, and even then the floor is $n(\log\log n)^2/\log n$. Whether that floor can actually be reached is **open**. Theorem C is the only obstacle we know, and we have no matching construction.

---

## 8. Escapes examined

**(a) Non-linear summaries — no gain [R, with a caveat].** Suppose each group $h$ publishes an arbitrary summary $\sigma_h(x_h)$, and every foreign linear input $\ell(x_h)$ that some factor reads must be a function of $\sigma_h(x_h)$. Then $\sigma_h$ refines the level sets of the linear map $L_h$ that stacks those forms. So $|\mathrm{im}\,\sigma_h|\ge|L_h(\{0,1\}^{S_h})|\ge2^{\operatorname{rank}L_h}$ by Lemma 1.1, and $L_h$ itself is the coarsest valid summary. *Caveat:* this assumes the summary must determine the factor **inputs**. A summary that determines only the factor **values** (or only their product) is not covered.

**(b) Sum-product over measurements — constant factor at most [S].** Instead of enumerating $\vec y$, one could eliminate the variables $y_h$ along a tree decomposition of the group-interaction graph. Groups $g$ and $h$ are adjacent when a clause of $g$ reads $y_h$. The cost becomes $(\#\text{values})^{O(\mathrm{tw})}$. Up to the zero-coefficient exceptions bounded in Theorem A, every cut edge of $G_F$ produces such a read. So a balanced separator of $w$ groups lifts to a separator of $G_F$ with at most $ws$ vertices. Expansion forces $\Omega(\varphi n/B)$ separator vertices, so $\mathrm{tw}=\Omega(\varphi n/(Bs))$, a constant fraction of the touched groups. Elimination therefore saves at most a constant factor in the exponent. We have not written out the full proof, including the per-group value counts.

**(c) Alternative recovery — Hermite counting [R for degree-only recovery].** Let $\mathcal P_d$ be the polynomials of degree at most $d$. Order-$J$ Taylor data at $t=1,\dots,Q$ determines $H(0)$ for every $H\in\mathcal P_d$ if and only if $QJ>d$. For "only if": if $QJ\le d$, then $W=\prod_t(T-t)^J\in\mathcal P_d$ has zero data and $W(0)\ne0$. So any recovery that uses only the degree bound $3nR$ needs $q\ge\lfloor3nR/Q\rfloor$, as in Section 2.5. What this does not exclude: recovery that exploits the *true* degree of $H$, or structure beyond degree, or that obtains order-$j$ Taylor coefficients without Lemma 3.3's $|E|\le j$ inclusion–exclusion. Appendix D's tapered Hermite recovery (Lemma D.5, with 2-adic valuations) is one example.

**(d) The remaining open door.** Every bound above uses the fact that the separation identity holds **for each assignment $x$**: every factor is computable from its own group plus the measurements, as in Definition 3.1(4) and Lemma 3.2. Identities that hold **only after summing over all assignments** — cancellations in $\sum_x$ such as Fourier or 2-adic cancellation, or tensor-network rewiring as in Appendix D — are not constrained by Theorems A–D. That is where any further improvement would have to come from. A second gap: Definition 3.1 designs whose measured subspaces vary with $t$ are covered only by Theorem A.

---

## 9. Experiments [E]

Everything in this section is computational evidence. None of it is proved. It comes from the scripts in `experiments/codes/`.

- **Ansatz.** The searches use a symmetric ansatz, not general $M_0$. Coordinates carry home vectors from a binary code $\mathcal C\subseteq\mathbb F_2^Q$ of dimension $k$. $M_0=\sum_S g_ST_S$ is translation-invariant, with $|S|\le R-1$, and the $g_S$ are found by a mixed-integer linear program (MILP) that minimises the number of exceptional characters $\tau$. Rows marked "not optimal" hit a time limit.
- **The cube reproduces the paper.** For $\mathcal C=\mathbb F_2^Q$ (the paper's cube), the MILP optimum matches the binomial-window ("interval") $\tau$ of Section 2.7. We re-ran this for $Q\le6$ and $d=R-1\le3$; all values agree, with two $Q=6$ cases solved to the same value but not proven optimal within the time limit. The script's notes claim agreement up to $Q=8$.
- **Other code families give constant factors that fade.** Families tried: cut codes of $K_{a,b}$, products $C_4\times C_4$ and $C_4\times C_6$, blow-ups $C_6[3]$, $C_6[4]$, $C_8[2]$, $C_8[3]$, Reed–Muller codes $\mathrm{RM}(2,5)$ and $\mathrm{RM}(2,6)$, lines of $AG(2,8)$, and 201 small random sparse parity-check codes. We compare the achieved $Q/R$ at $\tau\le0.05$ (and at $0.01$) with the cube's $Q/R$ at the same dimension $k$:

  | Code | $k$ | ratio at $\tau\le0.05$ | ratio at $\tau\le0.01$ |
  |---|---|---|---|
  | $K_{4,4}$ | 7 | 1.52 | 1.23 |
  | $K_{3,3}$ | 5 | 1.29 | — |
  | $K_{6,6}$ | 11 | 1.21 | 1.18 |
  | $K_{7,7}$ | 13 | 1.16 | 1.26 |
  | $K_{10,10}$ | 19 | 0.93 | 1.00 |
  | $K_{12,12}$ | 23 | 0.86 | 0.84 |
  | $\mathrm{RM}(2,5)$ | 16 | 0.94 | — |
  | $\mathrm{RM}(2,6)$ | 22 | 0.88 | 1.15 |

  The cycle products, the blow-ups and $AG(2,8)$ did not reach $\tau\le0.05$ within the values of $R$ tried. None of the random sparse codes beat the cube baseline at the same $k$. The gains are at most constant factors and are gone by about $k\approx19$.
- **Why this is expected.** Two reductions, used in the scripts and argued only informally there, explain the pattern:
  1. *Pull-back:* a code design pulls back to the full cube $\mathbb F_2^Q$, so $\tau_{\mathcal C}(d)\ge\tau^*_{\text{cube}}(Q,d)$.
  2. *Information set:* restricting to an information set $I$ gives $\tau_{\mathcal C}(d)\ge\tau^*_{\text{cube}}(k,d_I(d))$, where $d_I(d)=\max_{\mathrm{wt}(c)\le d}\mathrm{wt}_I(c)$.

  Both rely on the cube optimum being the interval design, which has been verified only for small $k$. They are consistent with Theorems C and D: a code changes the blocked fraction $\mu$ and the point count $|X|$ only by constants, and those enter the exponent only as constants.

---

## 10. What is proved and what was corrected

**Rigorous.** Lemmas 1.1, 1.2, 2.1, 4.1, 4.2, 5.1, 5.2 and 6.1; Theorems A–D; Corollaries E and F; escape (a) under its stated assumption; escape (c) for degree-only recovery. All constants ($1/24$, $3/64$, $e^{-16}/\sqrt2$, $8$, $16$, …) are unoptimised.

**Sketch.** Escape (b).

**Experimental.** All of Section 9, including both reductions mentioned there.

**Corrections to the brief that commissioned this note.**
1. Theorem C needs the hypothesis $8(R-1)\le\mu Q$; without it the rank can be 0.
2. The constant $\varphi n/32$ applies to the literal row count. For realised values it is $3\varphi n/(64B)$.
3. Theorem D needs $R\le\mu Q/4$ and $\mu Q\ge32$ to get an absolute constant $c$.
4. Theorem B's $(c/\tau)^{Q/(2R)}$ holds only for $Q\ge2R(R-1)$, and only for fully (or $K$-nearly) blocked coordinates.
5. Theorem A's threshold is $\varepsilon<\varphi/12$ in this note's normalisation; the brief had $\alpha/3$.
6. The corollary is in fact $\Omega(n(\log\log n)^2/\log n)$, slightly stronger than the stated $\Omega(n\log\log n/\log n)$. With $Q=O(\log n)$ it is $\Omega(n/\sqrt{\log n})$.
7. Theorems B–D do not cover Definition 3.1 designs with $t$-dependent measured subspaces.

## References

- OpenAI, *The Exponential Time Hypothesis is False*, 17 Oct 2026, revised 25 Oct 2026. Cited: Definition 2.1; Sections 2.4–2.7, (2.14)–(2.16); Definition 3.1; Lemmas 3.2–3.3; Proposition 3.4; Appendix D (Theorem D.1, Lemma D.5).
- G. L. Miller, R. Peng, S. C. Xu, *Parallel graph decompositions using random shifts*, SPAA 2013.
- B. Bollobás, *Random Graphs*, 2nd ed., Cambridge University Press, 2001 (expansion of random regular graphs).
