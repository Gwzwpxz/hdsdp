<TeXmacs|2.1.1>

<style|<tuple|generic|triangle-list>>

<\body>
  <\hide-preamble>
    \;

    <assign|a|<macro|<math-bf|a>>>

    <assign|b|<macro|<math-bf|b>>>

    <assign|c|<macro|<math-bf|c>>>

    <assign|d|<macro|<math-bf|d>>>

    <assign|e|<macro|<math-bf|e>>>

    <assign|1|<macro|<math-bf|1>>>

    <assign|A|<macro|<math-bf|A>>>

    <assign|x|<macro|<math-bf|x>>>

    <assign|y|<macro|<math-bf|y>>>

    <assign|s|<macro|<math-bf|s>>>

    <assign|u|<macro|<math-bf|u>>>

    <assign|0|<macro|<math-bf|0>>>

    <assign|dl|<macro|<math|<wide*|d|\<bar\>>>>>

    <assign|du|<macro|<math|<wide|d|\<bar\>>>>>

    <assign|w|<macro|<math-bf|w>>>
  </hide-preamble>

  <doc-data|<doc-title|Learning to Solve Linear
  Programming>|<doc-author|<author-data|<\author-affiliation>
    \;

    Wenzhi Gao

    \;
  </author-affiliation>|<\author-affiliation>
    <date|>
  </author-affiliation>>>>

  In this note, we discuss the potential to apply online algorithm to offline
  linear programming problems taking form of resource allocation. The main
  algorithm goes one or multiple rounds through the LP data and outputs both
  primal and dual estimates for the original LP.\ 

  Our theoretical results utilizes the error bound condition of optimization
  over bounded convex polytope and the contribution is two fold: 1) for
  online LP, we leverage a stability analysis from learning theory to derive
  an <math|\<cal-O\><around*|(|<wide|\<mu\>|\<bar\>><sqrt|<frac|m|n>>|)>>
  convergence result of dual convergence <math|<around*|\<\|\|\>|<y><rsub|n>-<y><rsup|\<ast\>>|\<\|\|\>>>,
  where <math|<wide|\<mu\>|\<bar\>>> is the error bound parameter associated
  with problem data distribution. 2) for online/offline LP, using the error
  bound condition, we show that an adaptive subgradient-based online
  algorithm equipped with a restart mechamism achieves
  <math|\<cal-O\><around*|(|m log n|)>> regret. Both results match the
  current lower bound in online LP setting up to a constant factor and do not
  explicitly assume the non-generacy condition. Practically we show that the
  subgradient-based online algorithm can be efficiently used to approximately
  solve real-life LPs. Especially, our empirical results demonstrate that our
  online algorithm can be efficiently applied in a sifting (column
  generation) procedure for large-scale linear programs.

  <with|font-series|bold|TODO>: Exploit the chance of stochastic proximal
  point/Bregman divergence.

  <section|Preliminaries and backgrounds>

  We are interested in primal-dual LPs given by

  <\eqnarray*>
    <tformat|<table|<row|<cell|max<rsub|<x>>>|<cell|<c><rsup|\<top\>><x>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<A><x>\<leq\><b>>|<cell|>>|<row|<cell|>|<cell|<0>\<leq\><x>\<leq\><u>>|<cell|>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<y>>>|<cell|<b><rsup|\<top\>><y>+<u><rsup|\<top\>><s>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|<A><rsup|\<top\>><y>+<s>\<geq\><c>>|<cell|>>|<row|<cell|>|<cell|<y>,<s>\<geq\>0.>|<cell|>>>>
  </eqnarray*>

  and our method is based on the following finite-sum reformulation

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<y>\<geq\><0>>>|<cell|<b><rsup|\<top\>><y>+<big|sum><rsub|i=1><rsup|n>u<rsub|i><around*|[|<c><rsub|i>-<a><rsub|i><rsup|\<top\>><y>|]><rsub|+>>|<cell|=:<big|sum><rsub|i=1><rsup|n>f<around*|(|<y>,\<xi\><rsub|i>|)>,>>>>
  </eqnarray*>

  where <math|\<xi\><rsub|i>=<around*|(|c<rsub|i>,<a><rsub|i>|)>> and
  <math|f<around*|(|<y>,\<xi\><rsub|i>|)>\<assign\><d><rsup|\<top\>><y>+u<rsub|i><around*|[|c<rsub|i>-<a><rsub|i><rsup|\<top\>><y>|]><rsub|+>,<d>=<b>/n>
  and without loss of geneality we take <math|<u>=<1>>. Througout the
  manuscript we make the following statistical assumptions

  <\render-theorem|Assumption 1>
    <em|(Stochastic input model)>

    <\em>
      <\itemize>
        <item>(i.i.d. sampling).\ 

        LP columns are drawn i.i.d. from certain distribution
        <math|\<cal-D\>>

        <item>(Bounded coefficient).

        There exists constant <math|<wide|c|\<bar\>>> and
        <math|<wide|a|\<bar\>>> such that
        <math|<around*|\<\|\|\>|<c>|\<\|\|\>><rsub|\<infty\>>\<leq\><wide|c|\<bar\>>>
        and <math|<around*|\<\|\|\>|<a><rsub|j>|\<\|\|\>>\<leq\><wide|a|\<bar\>>,\<forall\>j\<in\><around*|[|n|]>>.

        <item>(Linear growth of RHS).

        There exist constants <math|<wide|d|\<bar\>>,<wide*|d|\<bar\>>> such
        that <math|<wide*|d|\<bar\>><e>\<leq\><d>\<leq\><wide|d|\<bar\>><e>>.
      </itemize>
    </em>
  </render-theorem>

  The assumptions above apply to both online and offline LPs but we note that
  they are significantly different since

  <\enumerate>
    <item>In the offline LP setting, <math|\<cal-D\>> is taken to be a
    discrete uniform distribution over the <math|n> LP columns while in the
    online setting, <math|\<cal-D\>> could be a multivariate continuous
    distribution.

    <item>In the offline LP setting, the bounded coefficients are trivially
    accessible scanning through the LP coefficients and it's possible to
    design algorithms based on these bounds. However, in the online setting,
    although we assume the existence of the bounds, we do not know them
    before hand. However this difference does not matter as we employ a
    parameter-free adaptive framework.
  </enumerate>

  <subsection|Related works>

  This work is motivated by <cite|10.48550/arxiv.1909.05499> and
  <cite|10.48550/arxiv.2003.02513>, where the references therein have
  introduced the trend of research before it. Therefore we focus on the
  further researches on online LPs, learning theory and first order methods
  for linear programs.\ 

  <with|font-series|bold|Dual convergence>

  First, <cite|10.48550/arxiv.1909.05499> discusses dual convergence of
  online LP and obtains <math|\<cal-O\><around*|(|<sqrt|<frac|log log
  n|n>>|)>> bound under non-degeneracy assumption.
  <cite|10.2139/ssrn.3497056> proves that if <math|<a>\<geq\><0>>, then
  <math|<sqrt|n><around*|(|<y><rsub|n><rsup|\<ast\>>-<y><rsup|\<ast\>>|)>>
  converges in distribution to a normal random vector.
  <cite|10.48550/arxiv.2210.07996> imposes growth condition on
  <math|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y>,\<xi\>|)>|]>>
  directly and relaxes the non-degeneracy assumption.

  <with|font-series|bold|Online Regret bound>

  Since <cite|10.48550/arxiv.1909.05499> and
  <cite|10.48550/arxiv.2003.02513>, a line of research has improved over the
  regret bound and dependency over degeneracy. Specially
  <cite|10.48550/arxiv.2011.10124> considers a general online mirror descent
  framework, getting better dependency over <math|m>;
  <cite|10.2139/ssrn.3497056> improves the <math|\<cal-O\><around*|(|log log
  n log n|)>> rate to <math|\<cal-O\><around*|(|log n|)>> due to his shaper
  analysis of dual convergence; he also provides a lower bound showing that
  <math|\<cal-O\><around*|(|log n|)>> is optimal for the dual-based method.
  Similarly <cite|10.48550/arxiv.2210.07996> has improved dual convergence
  and gets <math|\<cal-O\><around*|(|log n|)>> regret without strictly
  satisfying the non-degeneracy condition.

  <big-table|<block|<tformat|<cwith|1|-1|1|-1|cell-halign|c>|<cwith|1|-1|1|-1|cell-tborder|1ln>|<cwith|1|-1|1|-1|cell-bborder|1ln>|<cwith|1|-1|1|-1|cell-lborder|0ln>|<cwith|1|-1|1|-1|cell-rborder|0ln>|<cwith|4|4|1|5|cell-halign|c>|<cwith|4|4|1|5|cell-tborder|1ln>|<cwith|4|4|1|5|cell-bborder|1ln>|<cwith|4|4|1|5|cell-lborder|0ln>|<cwith|4|4|1|5|cell-rborder|0ln>|<table|<row|<cell|Paper>|<cell|Briefing>|<cell|Regret>|<cell|Non-degeneracy>|<cell|Problem
  type>>|<row|<cell|<cite|10.48550/arxiv.1909.05499>>|<cell|LP-based>|<cell|<math|\<cal-O\><around*|(|log
  log n log n|)>>>|<cell|Yes>|<cell|LP online resource
  allocation>>|<row|<cell|<cite|10.48550/arxiv.2003.02513>>|<cell|Simple and
  fast>|<cell|<math|\<cal-O\><around*|(|<sqrt|n>|)>>>|<cell|No>|<cell|LP
  online resource allocation>>|<row|<cell|<cite|10.48550/arxiv.2011.10124>>|<cell|Mirror
  descent>|<cell|<math|\<cal-O\><around*|(|<sqrt|n>|)>>>|<cell|No>|<cell|General
  online resource allocation>>|<row|<cell|<cite|10.2139/ssrn.3497056>>|<cell|LP-based>|<cell|<math|\<cal-O\><around*|(|log
  n|)>>>|<cell|Yes>|<cell|LP online resource
  allocation>>|<row|<cell|<cite|10.48550/arxiv.2210.07996>>|<cell|LP-based>|<cell|<math|\<cal-O\><around*|(|log
  n|)>>>|<cell|Partly>|<cell|LP online resource allocation>>>>>|Research to
  now>

  <with|font-series|bold|Stability>

  <cite|hardt2016train> and <cite|shalev2010learnability> considers
  generalization, which is dual convergence in the online setting, using a
  neat stability argument.

  <with|font-series|bold|Error bound condition>

  Several first order methods have been analyzed under the second order
  growth condition, or more generally, the error bound condition
  <cite|10.48550/arxiv.1512.03107><cite|10.48550/arxiv.1805.04577>. We remark
  that in <cite|10.48550/arxiv.1805.04577>, an adaptive algorithm is proposed
  for finite-sum subgradient method on Lipschitz functions, improving
  convergence from <math|\<cal-O\><around*|(|<frac|1|<sqrt|n>>|)>> to
  <math|\<cal-O\><around*|(|<frac|log n|n>|)>>, which is the motivation of
  our idea.

  <\question>
    It seems that non-degeneracy condition is always used together with
    LP-based methods to reach the lower bound <math|\<cal-O\><around*|(|log
    n|)>>. What if <cite|10.48550/arxiv.2003.02513> considers the stochastic
    input model under the non-degenerate assumption (or, with the quadratic
    growth condition)? Will <math|\<cal-O\><around*|(|log n|)>> regret be
    possible, or what's the major difficulty (e.g., unknown growth
    parameter)?
  </question>

  <\question>
    Will a faster non-LP-based first-order dual method be benefitial to the
    simple and fast online algorithm? Or, does non LP-based method use dual
    convergence?
  </question>

  <\question>
    On the choice of potential function. Why use
    <math|<around*|\<\|\|\>|\<cdummy\>|\<\|\|\>><rsup|2>>?
  </question>

  <\question>
    Is it meaningful improving the <math|log n> term in random permutation?
  </question>

  <section|Error bound of dual LP>

  In this section, we focus on the error bound condition of dual LP. Namely
  we consider the following two problems

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<y>\<geq\><0>>>|<cell|<frac|1|n><big|sum><rsub|i=1><rsup|n>f<around*|(|<y>,\<xi\><rsub|i>|)>=<d><rsup|\<top\>><y>+<frac|1|n><big|sum><rsub|i=1><rsup|n><around*|[|<c><rsub|i>-<a><rsub|i><rsup|\<top\>><y>|]><rsub|+>>|<cell|<around*|(|\<cal-P\><rsub|n>|)>>>|<row|<cell|min<rsub|<y>\<geq\><0>>>|<cell|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y>,\<xi\>|)>|]>=<d><rsup|\<top\>><y>+\<bbb-E\><rsub|\<xi\>><around*|[|<around*|[|<c>-<a><rsup|\<top\>><y>|]><rsub|+>|]>>|<cell|<around*|(|\<cal-P\>|)>>>>>
  </eqnarray*>

  We need the following result from error bound literature (e.g.,
  <cite|10.48550/arxiv.1805.04577>).

  <\lemma>
    <\em>
      (Error bound of piecewise-linear function).\ 

      Convex optimization problem

      <\eqnarray*>
        <tformat|<table|<row|<cell|min<rsub|<y>>>|<cell|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y>,\<xi\>|)>|]>>|<cell|>>|<row|<cell|<text|subject
        to>>|<cell|<y>\<in\>\<cal-Y\>>|<cell|>>>>
      </eqnarray*>

      satisfies quadratic error bound condition if
      <math|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y>,\<xi\>|)>|]>> is
      piecewise-linear and <math|\<cal-Y\>> is a bounded polyhedron. i.e.,
      there exists <math|\<alpha\>\<gtr\>0> such that

      <\equation*>
        <around*|\<\|\|\>|<y>-Proj<rsub|\<cal-Y\><rsup|\<ast\>>><around*|(|<y>|)>|\<\|\|\>><rsup|2>\<leq\>\<alpha\><around*|{|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y>,\<xi\>|)>|]>-\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|Proj<rsub|\<cal-Y\><rsup|\<ast\>>><around*|(|<y>|)>,\<xi\>|)>|]>|}>,
      </equation*>

      where <math|\<cal-Y\><rsup|\<ast\>>> is the set of optimal solutions to
      the problem and <math|Proj<rsub|\<cal-Y\><rsup|\<ast\>>><around*|(|<y>|)>>
      denotes the orthogonal projection.
    </em>
  </lemma>

  Now we show that for each finite <math|n\<gtr\>0> and realization of the LP
  <math|<around*|(|<A>,<c>|)>>, the error bound holds.

  <\lemma>
    <\em>
      (Error bound of <math|><math|\<cal-P\><rsub|n>>)

      Assume that <strong|A1> holds. Then given <math|n\<gtr\>0> and LP
      realization <math|<around*|(|<A>,<c>|)>>, there exists some
      <math|\<alpha\><around*|(|<A>,<c>|)>\<gtr\>0> such that the error bound
      condition holds.\ 
    </em>
  </lemma>

  <\proof>
    Since <math|<0>> is a feasible solution to <math|\<cal-P\><rsub|n>>, we
    know that for any <math|<y><rsup|\<ast\>>\<in\>\<cal-Y\><rsup|\<ast\>>>

    <\equation*>
      <dl><around*|\<\|\|\>|<y><rsup|\<ast\>>|\<\|\|\>><rsub|1>\<leq\><d><rsup|\<top\>><y><rsup|\<ast\>>\<leq\><frac|1|n><big|sum><rsub|i=1><rsup|n>f<around*|(|<y><rsup|\<ast\>>,\<xi\><rsub|i>|)>\<leq\><frac|1|n><big|sum><rsub|i=1><rsup|n>f<around*|(|<0>,\<xi\><rsub|i>|)>\<leq\><wide|c|\<bar\>>
    </equation*>

    and <math|<y><rsup|\<ast\>>\<in\>\<Delta\><around*|(|<frac|<wide|c|\<bar\>>|<dl>>|)>>
    a simplex. Also we know that <math|<frac|1|n><big|sum><rsub|i=1><rsup|n>f<around*|(|<y>,\<xi\>|)>>
    is sum of convex piece-wise linear functions and applying <strong|Lemma
    1> completes the proof.
  </proof>

  <\remark>
    If distribution of <math|\<xi\>> has finite support, then
    <math|\<cal-P\>> also satisfies error bound condition.
  </remark>

  The error bound condition is the key to our main results. Now we define the
  notion error bound constant for online LP.

  <\definition>
    <\em>
      (Average error bound)

      The average error bound of online LP over distribution <math|\<cal-D\>>
      and decision horizon <math|n> is defined by

      <\equation*>
        \<mu\><rsub|n>\<assign\>\<bbb-E\><rsub|<around*|{|<around*|(|<a><rsub|j>,c<rsub|j>|)>\<sim\>\<cal-D\>,j=1\<ldots\>n|}>><around*|[|\<alpha\><around*|(|<A>,<c>|)>|]>.
      </equation*>
    </em>
  </definition>

  Our analysis makes the following assumption.

  <strong|Assumption 2.> (Uniform bounded error bound)

  The error bound satisfies <math|\<mu\><rsub|n>\<less\><wide|\<mu\>|\<bar\>>>
  for any <math|n\<gtr\>0>.

  <with|font-series|bold|TODO: justify it.>

  <section|Improving dual convergence>

  In this section, we attempt to improve the dual convergence result from
  online previous work.

  <with|font-series|bold|Use an algorithm that removes degeneracy.> e.g. add
  quadratic perturbation.

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<y>>>|<cell|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y>,\<xi\>|)>|]>+<frac|\<delta\>|2><around*|\<\|\|\>|<y>|\<\|\|\>><rsup|2>>|<cell|>>|<row|<cell|min<rsub|<y>>>|<cell|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y>,\<xi\>|)>|]>>|<cell|>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<y><rsub|1>>|<cell|=>|<cell|arg min<rsub|<y>>
    <around*|{|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y>,\<xi\>|)>|]>+<frac|\<delta\>|2><around*|\<\|\|\>|<y>|\<\|\|\>><rsup|2>|}>>>|<row|<cell|<y><rsub|2>>|<cell|\<in\>>|<cell|arg
    min<rsub|<y>> <around*|{|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y>,\<xi\>|)>|]>|}>>>|<row|<cell|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y><rsub|1>,\<xi\>|)>|]>>|<cell|\<leq\>>|<cell|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y><rsub|2>,\<xi\>|)>|]>+<frac|\<delta\>|2><around*|\<\|\|\>|<y><rsub|2>|\<\|\|\>><rsup|2>-<frac|\<delta\>|2><around*|\<\|\|\>|<y><rsub|1>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y><rsub|2>,\<xi\>|)>|]>+<frac|\<delta\>\<Delta\>|2m>>>>>
  </eqnarray*>

  Given two LPs that share <math|<around*|{|\<xi\><rsub|1>,\<ldots\>,\<xi\><rsub|j-1>,\<xi\><rsub|j+1>,\<ldots\>,\<xi\><rsub|n>|}>>
  and differ only in <math|\<xi\><rsub|j>,\<xi\><rsub|j><rprime|'>>, denote

  <\eqnarray*>
    <tformat|<table|<row|<cell|<y><rsup|\<ast\>><rsub|1>>|<cell|\<assign\>>|<cell|arg
    min<rsub|<y>\<in\>\<Delta\>> <around*|{|<frac|1|n><big|sum><rsub|i\<neq\>j>f<around*|(|<y>,\<xi\><rsub|i>|)>+<frac|1|n>f<around*|(|<y>,\<xi\><rsub|j>|)>|}>>>|<row|<cell|<y><rsup|\<ast\>><rsub|2>>|<cell|\<assign\>>|<cell|arg
    min<rsub|<y>\<in\>\<Delta\>> <around*|{|<frac|1|n><big|sum><rsub|i\<neq\>j>f<around*|(|<y>,\<xi\><rsub|i>|)>+<frac|1|n>f<around*|(|<y>,\<xi\><rprime|'><rsub|j>|)>|}>>>>>
  </eqnarray*>

  Due to the error bound condition, we know that

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<mu\><rsup|-1><around*|\<\|\|\>|P<rsub|\<cal-Y\><rsup|\<ast\>>><around*|[|<y><rsup|\<ast\>><rsub|2>|]>-<y><rsup|\<ast\>><rsub|2>|\<\|\|\>><rsup|2>>|<cell|\<leq\>>|<cell|<around*|{|<frac|1|n><big|sum><rsub|i\<neq\>j>f<around*|(|<y><rsup|\<ast\>><rsub|2>,\<xi\><rsub|i>|)>+<frac|1|n>f<around*|(|<y><rsup|\<ast\>><rsub|2>,\<xi\><rsub|j>|)>|}>-<around*|{|<frac|1|n><big|sum><rsub|i\<neq\>j>f<around*|(|P<around*|[|<y><rsup|\<ast\>><rsub|2>|]>,\<xi\><rsub|i>|)>+<frac|1|n>f<around*|(|P<around*|[|<y><rsup|\<ast\>><rsub|2>|]>,\<xi\><rsub|j>|)>|}>>>|<row|<cell|\<mu\><rsup|-1><around*|\<\|\|\>|<y><rsub|1><rsup|\<ast\>>-<y><rsup|\<ast\>><rsub|2>|\<\|\|\>><rsup|2>>|<cell|\<leq\>>|<cell|<around*|{|<frac|1|n><big|sum><rsub|i\<neq\>j>f<around*|(|<y><rsub|1><rsup|\<ast\>>,\<xi\><rsub|i>|)>+<frac|1|n>f<around*|(|<y><rsub|1><rsup|\<ast\>>,\<xi\><rprime|'><rsub|j>|)>|}>-<around*|{|<frac|1|n><big|sum><rsub|i\<neq\>j>f<around*|(|<y><rsup|\<ast\>><rsub|2>,\<xi\><rsub|i>|)>+<frac|1|n>f<around*|(|<y><rsup|\<ast\>><rsub|2>,\<xi\><rprime|'><rsub|j>|)>|}>.>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<mu\><rsup|-1><around*|\<\|\|\>|<y><rsub|1><rsup|\<ast\>>-<y><rsup|\<ast\>><rsub|2>|\<\|\|\>><rsup|2>>|<cell|\<leq\>>|<cell|<around*|{|<frac|1|n><big|sum><rsub|i\<neq\>j>f<around*|(|<y><rsup|\<ast\>><rsub|2>,\<xi\><rsub|i>|)>+<frac|1|n>f<around*|(|<y><rsup|\<ast\>><rsub|2>,\<xi\><rsub|j>|)>|}>-<around*|{|<frac|1|n><big|sum><rsub|i\<neq\>j>f<around*|(|<y><rsub|1><rsup|\<ast\>>,\<xi\><rsub|i>|)>+<frac|1|n>f<around*|(|<y><rsub|1><rsup|\<ast\>>,\<xi\><rsub|j>|)>|}>>>|<row|<cell|\<mu\><rsup|-1><around*|\<\|\|\>|<y><rsub|1><rsup|\<ast\>>-<y><rsup|\<ast\>><rsub|2>|\<\|\|\>><rsup|2>>|<cell|\<leq\>>|<cell|<around*|{|<frac|1|n><big|sum><rsub|i\<neq\>j>f<around*|(|<y><rsub|1><rsup|\<ast\>>,\<xi\><rsub|i>|)>+<frac|1|n>f<around*|(|<y><rsub|1><rsup|\<ast\>>,\<xi\><rprime|'><rsub|j>|)>|}>-<around*|{|<frac|1|n><big|sum><rsub|i\<neq\>j>f<around*|(|<y><rsup|\<ast\>><rsub|2>,\<xi\><rsub|i>|)>+<frac|1|n>f<around*|(|<y><rsup|\<ast\>><rsub|2>,\<xi\><rprime|'><rsub|j>|)>|}>.>>>>
  </eqnarray*>

  Taking summation over the two relations,

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|2\<mu\><rsup|-1><around*|\<\|\|\>|<y><rsub|1><rsup|\<ast\>>-<y><rsup|\<ast\>><rsub|2>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|1|n><around*|{|f<around*|(|<y><rsup|\<ast\>><rsub|2>,\<xi\><rsub|j>|)>-f<around*|(|<y><rsub|1><rsup|\<ast\>>,\<xi\><rsub|j>|)>+f<around*|(|<y><rsub|1><rsup|\<ast\>>,\<xi\><rprime|'><rsub|j>|)>-f<around*|(|<y><rsup|\<ast\>><rsub|2>,\<xi\><rprime|'><rsub|j>|)>|}>.>>>>
  </eqnarray*>

  Recall that

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|f<around*|(|<y><rsub|1><rsup|\<ast\>>,\<xi\><rsub|j>|)>-f<around*|(|<y><rsup|\<ast\>><rsub|2>,\<xi\><rsub|j>|)>>>|<row|<cell|>|<cell|=>|<cell|<d><rsup|\<top\>><y><rsub|1><rsup|\<ast\>>+<around*|[|c<rsub|i>-<a><rsub|i><rsup|\<top\>><y><rsub|1><rsup|\<ast\>>|]><rsub|+>-<d><rsup|\<top\>><y><rsub|2><rsup|\<ast\>>-<around*|[|c<rsub|i>-<a><rsub|i><rsup|\<top\>><y><rsub|2><rsup|\<ast\>>|]><rsub|+>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|(|<wide|a|\<bar\>>+<du>|)><sqrt|m><around*|\<\|\|\>|<y><rsub|1><rsup|\<ast\>>-<y><rsub|2><rsup|\<ast\>>|\<\|\|\>>,>>>>
  </eqnarray*>

  we know that

  <\equation*>
    2\<mu\><rsup|-1><around*|\<\|\|\>|<y><rsub|1><rsup|\<ast\>>-<y><rsup|\<ast\>><rsub|2>|\<\|\|\>><rsup|2>\<leq\><frac|<around*|(|<wide|a|\<bar\>>+<du>|)><sqrt|m>|n><around*|\<\|\|\>|<y><rsub|1><rsup|\<ast\>>-<y><rsub|2><rsup|\<ast\>>|\<\|\|\>>.
  </equation*>

  Dividing both sides by <math|<around*|\<\|\|\>|<y><rsub|1><rsup|\<ast\>>-<y><rsup|\<ast\>><rsub|2>|\<\|\|\>>>,
  we have

  <\equation*>
    <around*|\<\|\|\>|<y><rsub|1><rsup|\<ast\>>-<y><rsup|\<ast\>><rsub|2>|\<\|\|\>>\<leq\><frac|2\<mu\><around*|(|<wide|a|\<bar\>>+<du>|)><sqrt|m>|n>
  </equation*>

  and observing that

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|sup<rsub|\<xi\>\<sim\>\<cal-D\>>
    <around*|{|f<around*|(|<y><rsub|1><rsup|\<ast\>>,\<xi\>|)>-f<around*|(|<y><rsub|2><rsup|\<ast\>>,\<xi\>|)>|}>>>|<row|<cell|>|<cell|=>|<cell|sup<rsub|<around*|(|<a>,c|)>\<sim\>\<cal-D\>><around*|{|<d><rsup|\<top\>><y><rsub|1><rsup|\<ast\>>+<around*|[|c-<a><rsup|\<top\>><y><rsub|1><rsup|\<ast\>>|]><rsub|+>-<d><rsup|\<top\>><y><rsub|2><rsup|\<ast\>>-<around*|[|c-<a><rsup|\<top\>><y><rsub|2><rsup|\<ast\>>|]><rsub|+>|}>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|(|<wide|a|\<bar\>>+<du>|)><sqrt|m><around*|\<\|\|\>|<y><rsub|1><rsup|\<ast\>>-<y><rsub|2><rsup|\<ast\>>|\<\|\|\>>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|2\<mu\><around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>m|n>>>>>
  </eqnarray*>

  By the theory of algorithm stability, we immediately have

  <\equation*>
    \<bbb-E\><around*|[|<around*|\||<big|sum><rsub|i=1><rsup|n>f<around*|(|<y>,\<xi\><rsub|i>|)>-\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y>,\<xi\>|)>|]>|\|>|]>\<leq\><frac|2\<mu\><around*|(|<wide|a|\<bar\>>+<du>|)>m|n>.
  </equation*>

  and by the above Lipschitzian propertym we successively deduce that

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<around*|\<\|\|\>|<y><rsup|\<ast\>>-<y><rsub|n><rsup|\<ast\>>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|\<mu\><around*|[|\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y><rsub|n><rsup|\<ast\>>,\<xi\>|)>|]>-\<bbb-E\><rsub|\<xi\>><around*|[|f<around*|(|<y><rsup|\<ast\>>,\<xi\>|)>|]>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|2\<mu\><rsup|2><around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>m|n>>>>>
  </eqnarray*>

  and finally we recover

  <\equation*>
    \<bbb-E\><around*|[|<around*|\<\|\|\>|<y><rsup|\<ast\>>-<y><rsub|n><rsup|\<ast\>>|\<\|\|\>>|]>\<leq\><frac|<sqrt|2m>\<mu\><around*|(|<wide|a|\<bar\>>+<du>|)>|<sqrt|n>>,
  </equation*>

  which removes <math|log> term from the original bound.

  <section|Simple and fast implicit update for online LP>

  In this section, we discuss the application of implicit update, say

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|(|<y><rsup|t+1>,x<rsub|t>|)>>|<cell|=>|<cell|arg
    min<rsub|<y>\<geq\><0>> <around*|{|<d><rsup|\<top\>><y>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>|]><rsub|+>+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y>-<y><rsup|t>|\<\|\|\>><rsup|2>|}>,>>>>
  </eqnarray*>

  in the context of online linear programming. By the result from online LP,
  first we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<around*|[|<d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>|]>-<around*|[|<d><rsup|\<top\>><y><rsup|t+1>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>|]>>>|<row|<cell|>|<cell|=>|<cell|<around*|[|<d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>|]>-<around*|[|<d><rsup|\<top\>><y><rsup|t+1>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>|]>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n><around*|[|<d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>|]>-c<rsub|t>x<rsub|t>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|t=1><rsup|n><around*|[|<d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>|]>-c<rsub|t>x<rsub|t>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n><around*|[|<d><rsup|\<top\>><y><rsup|t>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t>|]><rsub|+>|]>-<around*|{|<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>+<around*|(|<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|)>x<rsub|t>|}>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n><around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>+2<sqrt|m><wide|a|\<bar\>><around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|t>|\<\|\|\>>>>>>
  </eqnarray*>

  \;

  Given a new column <math|<around*|(|c<rsub|t>,<a><rsub|t>|)>>, we
  re-introduce the variable <math|s> and consider\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<y>>>|<cell|<d><rsup|\<top\>><y>+s+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y>-<y><rsup|t>|\<\|\|\>><rsup|2>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|s\<geq\>c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>>|<cell|x<rsub|t>>>|<row|<cell|>|<cell|s\<geq\>0,<y>\<geq\><0>>|<cell|v<rsub|t>,>>>>
  </eqnarray*>

  whose Lagrangian function and KKT condition are respectively given by

  <\eqnarray*>
    <tformat|<table|<row|<cell|L<around*|(|<y>,x,<s>,x<rsub|t>,v<rsub|t>|)>>|<cell|=>|<cell|<d><rsup|\<top\>><y>+s+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y>-<y><rsup|t>|\<\|\|\>><rsup|2>+x<rsub|t><around*|(|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>-s|)>-v<rsub|t>s-<y><rsup|\<top\>><w>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|s>|<cell|\<geq\>>|<cell|<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>|]><rsub|+>>>|<row|<cell|v<rsub|t>,x<rsub|t>,<y>,<w>>|<cell|\<geq\>>|<cell|<0>>>|<row|<cell|<d>+<frac|<y>-<y><rsup|t>|\<gamma\><rsub|t>>-<a><rsub|t>x<rsub|t>-<w>>|<cell|=>|<cell|<0>>>|<row|<cell|0\<leq\>x<rsub|t>>|<cell|\<leq\>>|<cell|1>>|<row|<cell|<y><rsup|\<top\>><w>>|<cell|=>|<cell|0>>|<row|<cell|x<rsub|t><around*|(|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>-s|)>>|<cell|=>|<cell|0>>|<row|<cell|<around*|(|1-x<rsub|t>|)>s>|<cell|=>|<cell|0>>>>
  </eqnarray*>

  and we have, implicitly, that

  <\equation*>
    <y><rsup|t+1>=<y><rsup|t>-\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)>+\<gamma\><rsub|t><w>
  </equation*>

  for some <math|<w>\<geq\><0>>. Then we successively deduce that\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>>|<cell|=>|<cell|<around*|\<\|\|\>|<y><rsup|t>-\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)>+\<gamma\><rsub|t><w><rsup|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>-\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>-2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>+\<gamma\><rsub|t><rsup|2><around*|\<\|\|\>|<d>-<a><rsub|t>x<rsub|t>|\<\|\|\>><rsup|2>>>>>
  </eqnarray*>

  and re-arranging the terms,

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>+\<gamma\><rsub|t><rsup|2><around*|\<\|\|\>|<d>-<a><rsub|t>x<rsub|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>+\<gamma\><rsub|t><rsup|2><around*|[|m<dl><rsup|2>+<around*|\<\|\|\>|<a><rsub|t>|\<\|\|\>><rsup|2>-2<around*|\<langle\>|<d>,<a><rsub|t>x<rsub|t>|\<rangle\>>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>+2\<gamma\><rsub|t><rsup|2>m<dl><rsup|2>+2\<gamma\><rsub|t><rsup|2><around*|\<\|\|\>|<a><rsub|t>|\<\|\|\>><rsup|2>>>>>
  </eqnarray*>

  where the second inequality follows from

  <\eqnarray*>
    <tformat|<table|<row|<cell|2<around*|\<langle\>|<d>,<a><rsub|t>x<rsub|t>|\<rangle\>>>|<cell|\<leq\>>|<cell|2<around*|\<\|\|\>|<d>|\<\|\|\>><rsup|2>+2<around*|\<\|\|\>|<a><rsub|t>x<rsub|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|2<around*|\<\|\|\>|<d>|\<\|\|\>><rsup|2>+2<around*|\<\|\|\>|<a><rsub|t>|\<\|\|\>><rsup|2>.>>>>
  </eqnarray*>

  Taking <math|\<gamma\><rsub|k>\<equiv\>\<gamma\>>, we telescope from
  <math|t=1,\<ldots\>,n> and\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><rsub|t=1><rsup|n><around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsub|t>>|<cell|\<leq\>>|<cell|\<gamma\><rsup|-1><big|sum><rsub|t=1><rsup|n><around*|[|<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>+2\<gamma\><rsup|2>m<dl><rsup|2>+2\<gamma\><rsup|2><around*|\<\|\|\>|<a><rsub|t>|\<\|\|\>><rsup|2>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|2\<gamma\><big|sum><rsub|t=1><rsup|n><around*|(|m<dl><rsup|2>+<big|sum><rsub|t=1><rsup|n><around*|\<\|\|\>|<a><rsub|t>|\<\|\|\>><rsup|2>|)>+\<gamma\><rsup|-1><around*|\<\|\|\>|<y><rsup|1>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|2\<gamma\><big|sum><rsub|t=1><rsup|n><around*|(|m<dl><rsup|2>+<big|sum><rsub|t=1><rsup|n><around*|\<\|\|\>|<a><rsub|t>|\<\|\|\>><rsup|2>|)>.>>>>
  </eqnarray*>

  Taking expectation on both sides gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><around*|[|<around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>|]>>>>>
  </eqnarray*>

  \;

  First we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|<d><rsup|\<top\>><y><rsup|t+1>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|t>|\<\|\|\>><rsup|2>>|<cell|\<leq\>>|<cell|<d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>-<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|\<ast\>>|\<\|\|\>><rsup|2>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<d><rsup|\<top\>><y><rsup|t+1>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|t>|\<\|\|\>><rsup|2>>|<cell|\<leq\>>|<cell|c<rsub|t>+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>-<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>>>|<row|<cell|<d><rsup|\<top\>><y><rsup|t+1>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|t>|\<\|\|\>><rsup|2>>|<cell|\<leq\>>|<cell|<d><rsup|\<top\>><y><rsup|t>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t>|]><rsub|+>-<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|t>|\<\|\|\>><rsup|2>>>>>
  </eqnarray*>

  \;

  \;

  <\bibliography|bib|tm-alpha|olp_ref>
    <\bib-list|15>
      <bibitem*|BLM20><label|bib-10.48550/arxiv.2011.10124>Santiago Balseiro,
      Haihao Lu<localize|, and >Vahab Mirrokni. <newblock>The Best of Many
      Worlds: Dual Mirror Descent for Online Allocation Problems.
      <newblock><with|font-shape|italic|ArXiv>, 2020.<newblock>

      <bibitem*|Bra19a><label|bib-10.2139/ssrn.3497056>Robert Bray.
      <newblock>Does the Multisecretary Problem Always Have Bounded Regret?
      <newblock><with|font-shape|italic|SSRN Electronic Journal>,
      2019.<newblock>

      <bibitem*|Bra19b><label|bib-bray2019logarithmic>Robert<nbsp>L Bray.
      <newblock>Logarithmic regret in multisecretary and online linear
      programming problems with continuous valuations.
      <newblock><with|font-shape|italic|ArXiv e-prints>, <localize|page >0,
      2019.<newblock>

      <bibitem*|CLY21a><label|bib-10.48550/arxiv.2101.11092>Guanting Chen,
      Xiaocheng Li<localize|, and >Yinyu Ye. <newblock>An Improved Analysis
      of LP-based Control for Revenue Management.
      <newblock><with|font-shape|italic|ArXiv>, 2021.<newblock>

      <bibitem*|CLY21b><label|bib-10.48550/arxiv.2110.14621>Guanting Chen,
      Xiaocheng Li<localize|, and >Yinyu Ye. <newblock>Fairer LP-based Online
      Allocation via Analytic Center. <newblock><with|font-shape|italic|ArXiv>,
      2021.<newblock>

      <bibitem*|HRS16><label|bib-hardt2016train>Moritz Hardt, Ben
      Recht<localize|, and >Yoram Singer. <newblock>Train faster, generalize
      better: stability of stochastic gradient descent.
      <newblock><localize|In ><with|font-shape|italic|International
      conference on machine learning>, <localize|pages >1225\U1234. PMLR,
      2016.<newblock>

      <bibitem*|JMZ22><label|bib-10.48550/arxiv.2210.07996>Jiashuo Jiang,
      Will Ma<localize|, and >Jiawei Zhang. <newblock>Degeneracy is OK:
      Logarithmic Regret for Network Revenue Management with Indiscrete
      Distributions. <newblock><with|font-shape|italic|ArXiv>,
      2022.<newblock>

      <bibitem*|JZ20><label|bib-10.48550/arxiv.2012.07933>Jiashuo
      Jiang<localize| and >Jiawei Zhang. <newblock>Online Resource Allocation
      with Stochastic Resource Consumption.
      <newblock><with|font-shape|italic|ArXiv>, 2020.<newblock>

      <bibitem*|LSY20><label|bib-10.48550/arxiv.2003.02513>Xiaocheng Li,
      Chunlin Sun<localize|, and >Yinyu Ye. <newblock>Simple and Fast
      Algorithm for Binary Integer and Online Linear Programming.
      <newblock><with|font-shape|italic|ArXiv>, 2020.<newblock>

      <bibitem*|LY19><label|bib-10.48550/arxiv.1909.05499>Xiaocheng
      Li<localize| and >Yinyu Ye. <newblock>Online Linear Programming: Dual
      Convergence, New Algorithms, and Regret Bounds.
      <newblock><with|font-shape|italic|ArXiv>, 2019.<newblock>

      <bibitem*|LZZ+18><label|bib-10.48550/arxiv.1805.04577>Mingrui Liu,
      Xiaoxuan Zhang, Lijun Zhang, Rong Jin<localize|, and >Tianbao Yang.
      <newblock>Fast Rates of ERM and Stochastic Approximation: Adaptive to
      Error Bound Conditions. <newblock><with|font-shape|italic|ArXiv>,
      2018.<newblock>

      <bibitem*|MMT+20><label|bib-10.1145/3357713.3384309>Konstantin
      Makarychev, Yury Makarychev, Madhur Tulsiani, Gautam Kamath, Julia
      Chuzhoy, Jan<nbsp>van<nbsp>den Brand, Yin<nbsp>Tat Lee, Aaron
      Sidford<localize|, and >Zhao Song. <newblock>Solving tall dense linear
      programs in nearly linear time. <newblock><with|font-shape|italic|Proceedings
      of the 52nd Annual ACM SIGACT Symposium on Theory of Computing>,
      <localize|pages >775\U788, 2020.<newblock>

      <bibitem*|SHW21><label|bib-10.48550/arxiv.2108.00331>Jinyan Su, Lijie
      Hu<localize|, and >Di Wang. <newblock>Faster Rates of Private
      Stochastic Convex Optimization. <newblock><with|font-shape|italic|ArXiv>,
      2021.<newblock>

      <bibitem*|SSSS10><label|bib-shalev2010learnability>Shai Shalev-Shwartz,
      Ohad Shamir, Nathan Srebro<localize|, and >Karthik Sridharan.
      <newblock>Learnability, stability and uniform convergence.
      <newblock><with|font-shape|italic|The Journal of Machine Learning
      Research>, 11:2635\U2670, 2010.<newblock>

      <bibitem*|YL15><label|bib-10.48550/arxiv.1512.03107>Tianbao
      Yang<localize| and >Qihang Lin. <newblock>RSG: Beating Subgradient
      Method without Smoothness and Strong Convexity.
      <newblock><with|font-shape|italic|ArXiv>, 2015.<newblock>
    </bib-list>
  </bibliography>

  <appendix|Improved Convergence of simple and fast online algorithm>

  In this section we duplicate the proof of simple and fast online algorithm
  and try to improve some of the algorithm constants. Recall that the
  algorithm updates by

  <\eqnarray*>
    <tformat|<table|<row|<cell|<y><rsup|t+1>>|<cell|=>|<cell|<around*|[|<y><rsup|t>+\<gamma\><rsub|t><around*|(|<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>-<d>|)><rsub|>|]><rsub|+>.>>>>
  </eqnarray*>

  and we have

  <\equation*>
    \<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>\<leq\><big|sum><rsub|t=1><rsup|n>\<bbb-E\><around*|[|<around*|(|<d>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsup|t>|}>|)><rsup|\<top\>><y><rsup|t>|]>.
  </equation*>

  The proof in <cite|10.48550/arxiv.2003.02513> adopts the potential function
  <math|\<phi\><around*|(|<y>|)>=<around*|\<\|\|\>|<y>|\<\|\|\>><rsup|2>> and
  deduces that

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<phi\><around*|(|<y><rsup|t+1>|)>-\<phi\><around*|(|<y><rsup|t>|)>>|<cell|=>|<cell|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<around*|\<\|\|\>|<around*|[|<y><rsup|t>+\<gamma\><rsub|t><around*|(|<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>-<d>|)><rsub|>|]><rsub|+>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>+\<gamma\><rsub|t><around*|(|<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>-<d>|)>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|-2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>|)><rsup|\<top\>><y><rsup|t>+\<gamma\><rsub|t><rsup|2><around*|\<\|\|\>|<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>-<d>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|-2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>|)><rsup|\<top\>><y><rsup|t>+m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>.>>>>
  </eqnarray*>

  Then telescoping over <math|t=1,\<ldots\>,n> gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><rsub|t=1><rsup|n>\<phi\><around*|(|<y><rsup|t+1>|)>-\<phi\><around*|(|<y><rsup|t>|)>>|<cell|=>|<cell|\<phi\><around*|(|<y><rsup|n+1>|)>-\<phi\><around*|(|<y><rsup|1>|)>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n><around*|[|-2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>|)><rsup|\<top\>><y><rsup|t>+m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>|]>>>>>
  </eqnarray*>

  and re-arranging the terms,

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<big|sum><rsub|t=1><rsup|n>2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>|)><rsup|\<top\>><y><rsup|t>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>-\<phi\><around*|(|<y><rsup|n+1>|)>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>>>>>
  </eqnarray*>

  and if we take <math|\<gamma\><rsub|t>\<equiv\>\<gamma\>>, then.

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><around*|[|<around*|(|<d><rsup|\<top\>>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsup|t>|}>|)><rsup|\<top\>><y><rsup|t>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\>n|2>.>>>>
  </eqnarray*>

  It can be readily seen that the bound of regret can be improved taking
  <math|\<gamma\>> arbitrarily small. But this is not feasible due to the
  huge constraint violation. Recall that to bound the constraint violation,
  we have

  <\equation*>
    <y><rsup|t+1>=<around*|[|<y><rsup|t>+\<gamma\><rsub|t><around*|(|<a>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>-<d>|)><rsub|>|]><rsub|+>\<geq\><y><rsup|t>+\<gamma\><around*|(|<a><rsub|t>x<rsub|t>-<d>|)>
  </equation*>

  and\ 

  <\equation*>
    <big|sum><rsub|t=1><rsup|n><a><rsub|t>x<rsub|t>\<leq\><b>+\<gamma\><rsup|-1><big|sum><rsub|t=1><rsup|n><around*|(|<y><rsup|t+1>-<y><rsup|t>|)>\<leq\><b>+\<gamma\><rsup|-1><y><rsup|t+1>
  </equation*>

  which finally gives

  <\equation*>
    \<bbb-E\><around*|[|<around*|\<\|\|\>|<around*|[|<A><x>-<b>|]><rsub|+>|\<\|\|\>>|]>\<leq\>\<gamma\><rsup|-1>\<bbb-E\><around*|[|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>>|]>
  </equation*>

  and we need to bound <math|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>>>. Next
  we look at the bound of <math|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>>>,
  where we successively deduce that\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>+\<gamma\><rsub|t><around*|(|<a><rsub|t>x<rsub|t>-<d>|)>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|-2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>+m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>+2\<gamma\><rsub|t><wide|c|\<bar\>>-2\<gamma\><rsub|t><dl><around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>>.>>>>
  </eqnarray*>

  On the other hand,\ 

  <\equation*>
    <around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>>\<leq\><around*|\<\|\|\>|<y><rsup|t>+\<gamma\><rsub|t><around*|(|<a><rsub|t>x<rsub|t>-<d>|)>|\<\|\|\>>\<leq\><around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>>+\<gamma\><rsub|t><around*|\<\|\|\>|<a><rsub|t>x<rsub|t>-<d>|\<\|\|\>>\<leq\><around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>>+\<gamma\><rsub|t><sqrt|m><around*|(|<wide|a|\<bar\>>+<du>|)>.
  </equation*>

  Here we need a threshold <math|\<tau\>> such that by induction
  <math|<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>>\<leq\>\<tau\>,\<forall\>t\<gtr\>0>.
  Since\ 

  <\equation*>
    <around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>>\<geq\><frac|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t>+2<wide|c|\<bar\>>|2<dl>>
  </equation*>

  would give <math|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>\<leq\><around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>>
  and <math|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>>\<leq\><around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>>+\<gamma\><rsub|t><sqrt|m><around*|(|<wide|a|\<bar\>>+<du>|)>>.\ 

  Given <math|\<gamma\><rsub|t>\<equiv\>\<gamma\>\<less\>1>, we have

  <\equation*>
    <around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>>\<leq\><frac|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\>|<dl>>+<frac|<wide|c|\<bar\>>|<dl>>+\<gamma\><sqrt|m><around*|(|<wide|a|\<bar\>>+<du>|)>
  </equation*>

  and\ 

  <\equation*>
    \<bbb-E\><around*|[|<around*|\<\|\|\>|<around*|[|<A><x>-<b>|]><rsub|+>|\<\|\|\>>|]>\<leq\>\<gamma\><rsup|-1>
    <around*|(|<frac|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\>|<dl>>+<frac|<wide|c|\<bar\>>|<dl>>+\<gamma\><sqrt|m><around*|(|<wide|a|\<bar\>>+<du>|)>|)>.
  </equation*>

  Now we take a look of the regret and constraint violation

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>>|<cell|\<leq\>>|<cell|<frac|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\>n|2>>>|<row|<cell|\<bbb-E\><around*|[|<around*|\<\|\|\>|<around*|[|<A><x>-<b>|]><rsub|+>|\<\|\|\>>|]>>|<cell|\<leq\>>|<cell|\<gamma\><rsup|-1>
    <around*|(|<frac|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\>|<dl>>+<frac|<wide|c|\<bar\>>|<dl>>+\<gamma\><sqrt|m><around*|(|<wide|a|\<bar\>>+<du>|)>|)>>>|<row|<cell|>|<cell|=>|<cell|<frac|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>|<dl>>+<sqrt|m><around*|(|<wide|a|\<bar\>>+<du>|)>+<frac|<wide|c|\<bar\>>|\<gamma\><dl>>>>>>
  </eqnarray*>

  If we wish we minimize the sum of <math|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>>
  and <math|\<bbb-E\><around*|[|<around*|\<\|\|\>|<around*|[|<A><x>-<b>|]><rsub|+>|\<\|\|\>>|]>>,

  <\equation*>
    \<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>+\<bbb-E\><around*|[|<around*|\<\|\|\>|<around*|[|<A><x>-<b>|]><rsub|+>|\<\|\|\>>|]>,
  </equation*>

  we have

  <\equation*>
    <frac|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\>n|2>+<frac|<wide|c|\<bar\>>|\<gamma\><dl>>\<leq\>2<around*|(|<frac|<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2><wide|c|\<bar\>>|2<dl>>|)><rsup|1/2>
  </equation*>

  taking

  <\equation*>
    \<gamma\>=<sqrt|<frac|2<wide|c|\<bar\>>|<dl><around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>
    m n>>\<sim\>\<cal-O\><around*|(|<frac|1|<sqrt|m n>>|)>.
  </equation*>

  Back to the convergence result, we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>>|<cell|\<leq\>>|<cell|2<around*|(|<frac|<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2><wide|c|\<bar\>>|2<dl>>|)><rsup|1/2><sqrt|m
    n>>>|<row|<cell|\<bbb-E\><around*|[|<around*|\<\|\|\>|<around*|[|<A><x>-<b>|]><rsub|+>|\<\|\|\>>|]>>|<cell|\<leq\>>|<cell|<frac|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>|<dl>>+<sqrt|m><around*|(|<wide|a|\<bar\>>+<du>|)>+<around*|(|<frac|<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2><wide|c|\<bar\>>|2<dl>>|)><rsup|1/2><sqrt|m
    n>,>>>>
  </eqnarray*>

  which reduces dependence over <math|m>.

  <with|font-series|bold|Use of new potential function>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<around*|[|<d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>|]>-<around*|[|<d><rsup|\<top\>><y><rsup|t>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t>|]><rsub|+>|]>>>|<row|<cell|>|<cell|=>|<cell|<around*|[|<d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>|]>-<around*|[|<d><rsup|\<top\>><y><rsup|t>+<around*|(|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t>|)>x<rsub|t>|]>>>|<row|<cell|>|<cell|=>|<cell|<around*|[|<d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>|]>-c<rsub|t>x<rsub|t>-<around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n><around*|[|<d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>|]>-c<rsub|t>x<rsub|t>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|t=1><rsup|n><around*|[|<d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>|]>-<around*|[|<d><rsup|\<top\>><y><rsup|t>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t>|]><rsub|+>|]>+<around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>>>>>
  </eqnarray*>

  Here we adopt the potential function <math|\<phi\><around*|(|<y>|)>\<assign\><around*|\<\|\|\>|<y>-<y><rsup|\<ast\>>|\<\|\|\>><rsup|2>>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|\<phi\><around*|(|<y><rsup|t+1>|)>-\<phi\><around*|(|<y><rsup|t>|)>>>|<row|<cell|>|<cell|=>|<cell|<around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|\<ast\>>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t>-<y><rsup|\<ast\>>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<around*|\<\|\|\>|<around*|[|<y><rsup|t>+\<gamma\><rsub|t><around*|(|<a><rsub|t>x<rsub|t>-<d>|)><rsub|>|]><rsub|+>-<y><rsup|\<ast\>>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t>-<y><rsup|\<ast\>>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>-<y><rsup|\<ast\>>-\<gamma\><rsub|t><around*|(|<a><rsub|t>x<rsub|t>-<d>|)>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t>-<y><rsup|\<ast\>>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|-2\<gamma\><rsub|t><around*|(|<a><rsub|t>x<rsub|t>-<d>|)><rsup|\<top\>><around*|(|<y><rsup|t>-<y><rsup|\<ast\>>|)>+\<gamma\><rsub|t><rsup|2><around*|\<\|\|\>|<a><rsub|t>x<rsub|t>-<d>|\<\|\|\>><rsup|2>>>>>
  </eqnarray*>

  \;

  \;

  Then telescoping over <math|t=1,\<ldots\>,n> gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><rsub|t=1><rsup|n>\<phi\><around*|(|<y><rsup|t+1>|)>-\<phi\><around*|(|<y><rsup|t>|)>>|<cell|=>|<cell|\<phi\><around*|(|<y><rsup|n+1>|)>-\<phi\><around*|(|<y><rsup|1>|)>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n><around*|[|-2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>|)><rsup|\<top\>><y><rsup|t>+m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>|]>>>>>
  </eqnarray*>

  and re-arranging the terms,

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<big|sum><rsub|t=1><rsup|n>2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>|)><rsup|\<top\>><y><rsup|t>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>-\<phi\><around*|(|<y><rsup|n+1>|)>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>>>>>
  </eqnarray*>

  and if we take <math|\<gamma\><rsub|t>\<equiv\>\<gamma\>>, then.

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><around*|[|<around*|(|<d><rsup|\<top\>>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsup|t>|}>|)><rsup|\<top\>><y><rsup|t>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\>n|2>.>>>>
  </eqnarray*>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|1.1|2>>
    <associate|auto-3|<tuple|1|2>>
    <associate|auto-4|<tuple|2|3>>
    <associate|auto-5|<tuple|3|4>>
    <associate|auto-6|<tuple|4|6>>
    <associate|auto-7|<tuple|4|7>>
    <associate|auto-8|<tuple|A|8>>
    <associate|bib-10.1145/3357713.3384309|<tuple|MMT+20|8>>
    <associate|bib-10.2139/ssrn.3497056|<tuple|Bra19a|7>>
    <associate|bib-10.48550/arxiv.1512.03107|<tuple|YL15|8>>
    <associate|bib-10.48550/arxiv.1805.04577|<tuple|LZZ+18|8>>
    <associate|bib-10.48550/arxiv.1909.05499|<tuple|LY19|8>>
    <associate|bib-10.48550/arxiv.2003.02513|<tuple|LSY20|8>>
    <associate|bib-10.48550/arxiv.2011.10124|<tuple|BLM20|7>>
    <associate|bib-10.48550/arxiv.2012.07933|<tuple|JZ20|8>>
    <associate|bib-10.48550/arxiv.2101.11092|<tuple|CLY21a|7>>
    <associate|bib-10.48550/arxiv.2108.00331|<tuple|SHW21|8>>
    <associate|bib-10.48550/arxiv.2110.14621|<tuple|CLY21b|8>>
    <associate|bib-10.48550/arxiv.2210.07996|<tuple|JMZ22|8>>
    <associate|bib-bray2019logarithmic|<tuple|Bra19b|7>>
    <associate|bib-hardt2016train|<tuple|HRS16|8>>
    <associate|bib-shalev2010learnability|<tuple|SSSS10|8>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      10.48550/arxiv.1909.05499

      10.48550/arxiv.2003.02513

      10.48550/arxiv.1909.05499

      10.2139/ssrn.3497056

      10.48550/arxiv.2210.07996

      10.48550/arxiv.1909.05499

      10.48550/arxiv.2003.02513

      10.48550/arxiv.2011.10124

      10.2139/ssrn.3497056

      10.48550/arxiv.2210.07996

      10.48550/arxiv.1909.05499

      10.48550/arxiv.2003.02513

      10.48550/arxiv.2011.10124

      10.2139/ssrn.3497056

      10.48550/arxiv.2210.07996

      hardt2016train

      shalev2010learnability

      10.48550/arxiv.1512.03107

      10.48550/arxiv.1805.04577

      10.48550/arxiv.1805.04577

      10.48550/arxiv.2003.02513

      10.48550/arxiv.1805.04577

      10.48550/arxiv.2003.02513
    </associate>
    <\associate|table>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Research to
      now>|<pageref|auto-3>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Preliminaries
      and backgrounds> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>Related works
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Error
      bound of dual LP> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Improving
      dual convergence> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Simple
      and fast implicit update for online LP>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Appendix
      A<space|2spc>Improved Convergence of simple and fast online algorithm>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>