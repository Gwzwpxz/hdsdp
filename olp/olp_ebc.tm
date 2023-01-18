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

  <doc-data|<doc-title|An improved analysis for dual subgradient online
  LP>|<doc-author|<author-data|<\author-affiliation>
    \;
  </author-affiliation>|<\author-affiliation>
    <date|>
  </author-affiliation>>>>

  <abstract-data|<\abstract>
    In this note we 1) give an improved analysis of simple and fast online
    algorithm which reduces both constraint violation and regret to
    <math|\<cal-O\><around*|(|<sqrt|m n>|)>>. 2) we try to develop a
    convergence analysis of simple and fast online algorithm under local
    error bound condition.
  </abstract>>

  In this note, we try to improve the analysis of dual subgradient-based
  online algorithm.

  <\eqnarray*>
    <tformat|<table|<row|<cell|x<rsub|t>>|<cell|=>|<cell|\<bbb-I\><around*|{|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>\<gtr\>0|}>>>|<row|<cell|<y><rsup|t+1>>|<cell|=>|<cell|<around*|[|<y><rsup|t>-\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)>|]><rsub|+><eq-number>>>>>
  </eqnarray*>

  The proof from [1] chooses suboptimal stepsize and achieves
  <math|\<cal-O\><around*|(|m<sqrt|n>|)>> violation.

  <with|font-series|bold|Notations>

  We use <math|\<bbb-E\><rsub|t><around*|[|\<cdummy\>|]>> to denote
  <math|\<bbb-E\><around*|[|\<cdummy\>\|\<sigma\><rsub|t>|]>>, where
  <math|\<sigma\><rsub|t>> is the <math|\<sigma\>>-algebra generating the
  past iterates. The coefficient bounds and assumptions are consistent with
  the original paper, namely <math|<wide|a|\<bar\>>,<dl>,<du>,<wide|c|\<bar\>>>.

  <section|Improved convergence of simple and fast online algorithm>

  In this section we duplicate the proof of simple and fast online algorithm
  [1] and try to improve some of the algorithm constants.\ 

  [1] Li, X., Sun, C. & Ye, Y. Simple and Fast Algorithm for Binary Integer
  and Online Linear Programming.<nbsp><with|font-shape|italic|Arxiv><nbsp>(2020)
  doi:10.48550/arxiv.2003.02513.

  <hrule>

  <with|font-series|bold|Regret bound>

  Recall that the algorithm updates by

  <\eqnarray*>
    <tformat|<table|<row|<cell|<y><rsup|t+1>>|<cell|=>|<cell|<around*|[|<y><rsup|t>+\<gamma\><rsub|t><around*|(|<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>-<d>|)><rsub|>|]><rsub|+>.>>>>
  </eqnarray*>

  and we have

  <\equation*>
    \<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>\<leq\><big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<around*|(|<d>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsup|t>|}>|)><rsup|\<top\>><y><rsup|t>|]>.
  </equation*>

  The original proof [1] adopts the potential function
  <math|\<phi\><around*|(|<y>|)>=<around*|\<\|\|\>|<y>|\<\|\|\>><rsup|2>> and
  deduces that

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<phi\><around*|(|<y><rsup|t+1>|)>-\<phi\><around*|(|<y><rsup|t>|)>>|<cell|=>|<cell|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<around*|\<\|\|\>|<around*|[|<y><rsup|t>+\<gamma\><rsub|t><around*|(|<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>-<d>|)><rsub|>|]><rsub|+>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>+\<gamma\><rsub|t><around*|(|<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>-<d>|)>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|-2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>|)><rsup|\<top\>><y><rsup|t>+\<gamma\><rsub|t><rsup|2><around*|\<\|\|\>|<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>-<d>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|-2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>|)><rsup|\<top\>><y><rsup|t>+m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>.>>>>
  </eqnarray*>

  Then we take expectation and telescope over <math|t=1,\<ldots\>,n> to
  obtain

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|\<phi\><around*|(|<y><rsup|t+1>|)>-\<phi\><around*|(|<y><rsup|t>|)>|]>>|<cell|=>|<cell|\<bbb-E\><around*|[|\<phi\><around*|(|<y><rsup|n+1>|)>|]>-\<phi\><around*|(|<y><rsup|1>|)>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|-2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>|)><rsup|\<top\>><y><rsup|t>+m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>|]>>>>>
  </eqnarray*>

  and re-arranging the terms,

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<big|sum><rsub|t=1><rsup|n>2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsub|t>|}>|)><rsup|\<top\>><y><rsup|t>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>-\<bbb-E\><rsub|t><around*|[|\<phi\><around*|(|<y><rsup|n+1>|)>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>>>>>
  </eqnarray*>

  and if we take <math|\<gamma\><rsub|t>\<equiv\>\<gamma\>>, then.

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<around*|(|<d><rsup|\<top\>>-<a><rsub|t>\<bbb-I\><around*|{|c<rsub|t>\<gtr\><a><rsub|t><rsup|\<top\>><y><rsup|t>|}>|)><rsup|\<top\>><y><rsup|t>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\>n|2>.>>>>
  </eqnarray*>

  It can be readily seen that the bound of regret can be improved taking
  <math|\<gamma\>> arbitrarily small. But this is not feasible due to the
  huge constraint violation.\ 

  <hrule>

  <with|font-series|bold|Constraint violation>

  Recall that to bound the constraint violation, we have

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
    <tformat|<table|<row|<cell|>|<cell|>|<cell|\<bbb-E\><rsub|t><around*|[|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>|]>-<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|\<bbb-E\><rsub|t><around*|[|<around*|\<\|\|\>|<y><rsup|t>+\<gamma\><rsub|t><around*|(|<a><rsub|t>x<rsub|t>-<d>|)>|\<\|\|\>><rsup|2>|]>-<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|\<bbb-E\><rsub|t><around*|[|-2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>|]>+m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsub|t><rsup|2>+2\<gamma\><rsub|t><wide|c|\<bar\>>-2\<gamma\><rsub|t><dl><around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>>.>>>>
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

  <hrule>

  <with|font-series|bold|Balancing regret and violation>

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
    <frac|m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\>n|2>+<frac|<wide|c|\<bar\>>|\<gamma\><dl>>\<leq\>2<around*|(|<frac|<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2><wide|c|\<bar\>>|2<dl>>|)><rsup|1/2><sqrt|m
    n>
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

  <section|Incorporating local second order growth>

  Let me informally assume that

  <\equation*>
    <big|sum><rsub|t=1><rsup|n><d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>-<big|sum><rsub|t=1><rsup|n><d><rsup|\<top\>><y>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>|]><rsub|+>\<leq\>-<frac|\<mu\>|2><around*|\<\|\|\>|<y>-<y><rsup|\<ast\>>|\<\|\|\>><rsup|2>
  </equation*>

  Since the local growth condition is relevant to set of optimal solutions.\ 

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

  <\equation*>
    <big|sum><rsub|t=1><rsup|n><d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>\<leq\><big|sum><rsub|t=1><rsup|n><d><rsup|\<top\>><y><rsup|t>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>|]><rsub|+>-<frac|\<mu\>|2><around*|\<\|\|\>|<y>-<y><rsup|\<ast\>>|\<\|\|\>><rsup|2>
  </equation*>

  \;

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>|]>>|<cell|\<leq\>>|<cell|n
    f<around*|(|<y><rsup|\<ast\>>|)>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|t=1><rsup|n>
    f<around*|(|<y><rsup|\<ast\>>|)>-<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t>
    <around*|[|f<rsub|n><around*|(|<y><rsup|t>|)>|]>+<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t>
    <around*|[|f<rsub|n><around*|(|<y><rsup|t>|)>|]>>>|<row|<cell|>|<cell|=>|<cell|<around*|{|<big|sum><rsub|t=1><rsup|n>
    f<around*|(|<y><rsup|\<ast\>>|)>-\<bbb-E\><rsub|t><around*|[|f<around*|(|<y><rsup|t>|)>|]>|}>+<around*|{|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|f<around*|(|<y><rsup|t>|)>|]>-<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|f<rsub|t><around*|(|<y><rsup|t>|)>|]>|}>+<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|f<rsub|t><around*|(|<y><rsup|t>|)>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|-<frac|\<mu\>|2><around*|\<\|\|\>|<y><rsup|t>-<y><rsup|\<ast\>>|\<\|\|\>><rsup|2>>>>>
  </eqnarray*>

  \;
</body>

<\initial>
  <\collection>
    <associate|font-base-size|10>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|2|3>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Improved
      convergence of simple and fast online algorithm>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Incorporating
      local second order growth> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>