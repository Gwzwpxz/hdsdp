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

  <doc-data|<doc-title|Implicit online LP>|<doc-author|<author-data|<\author-affiliation>
    \;
  </author-affiliation>|<\author-affiliation>
    <date|>
  </author-affiliation>>>>

  <abstract-data|<abstract|In this note we try to develop an online algorithm
  using implicit proximal point updates. The regret and constraint violation
  of the method are both <math|\<cal-O\><around*|(|<sqrt|m n>|)>> with
  potential better performance if maximum-norm dual optimal solution has
  relatively small size.>>

  In this note, we discuss an implicit dual update model for online LP with
  stochastic input, say

  <\eqnarray*>
    <tformat|<table|<row|<cell|<y><rsup|t+1>>|<cell|=>|<cell|arg
    min<rsub|<y>\<geq\><0>> <around*|{|<d><rsup|\<top\>><y>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>|]><rsub|+>+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y>-<y><rsup|t>|\<\|\|\>><rsup|2>|}><eq-number>>>>>
  </eqnarray*>

  The intuition of the algorithm is that <math|<d><rsup|\<top\>><y>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>|]><rsub|+>>
  would preserve more information than linearizing it using subgradient and
  thus is more likely to behave better in practice. This updating formula
  corresponds to stochastic proximal point, or implicit online algorithm in
  literature.\ 

  <with|font-series|bold|Notations>

  We use <math|\<bbb-E\><rsub|t><around*|[|\<cdummy\>|]>> to denote
  <math|\<bbb-E\><around*|[|\<cdummy\>\|\<sigma\><rsub|t>|]>>, where
  <math|\<sigma\><rsub|t>> is the <math|\<sigma\>>-algebra generating the
  past iterates. The coefficient bounds and assumptions are consistent with
  the original paper, namely <math|<wide|a|\<bar\>>,<dl>,<du>,<wide|c|\<bar\>>>.

  <hrule>

  <with|font-series|bold|Properties of implicit update>

  To derive the results we need, we first analyze the optimization problem
  (1) and derive basic properties of the implicit update. Given a new column
  <math|<around*|(|c<rsub|t>,<a><rsub|t>|)>>, we re-introduce the variable
  <math|s> and consider\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<y>>>|<cell|<d><rsup|\<top\>><y>+s+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y>-<y><rsup|t>|\<\|\|\>><rsup|2>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|s\<geq\>c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>>|<cell|x<rsub|t>>>|<row|<cell|>|<cell|s\<geq\>0,<y>\<geq\><0>>|<cell|v<rsub|t>,<w>>>>>
  </eqnarray*>

  whose Lagrangian function and KKT condition are respectively given by

  <\eqnarray*>
    <tformat|<table|<row|<cell|L<around*|(|<y>,x,<s>,x<rsub|t>,v<rsub|t>|)>>|<cell|=>|<cell|<d><rsup|\<top\>><y>+s+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y>-<y><rsup|t>|\<\|\|\>><rsup|2>+x<rsub|t><around*|(|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>-s|)>-v<rsub|t>s-<y><rsup|\<top\>><w>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|s>|<cell|\<geq\>>|<cell|<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>|]><rsub|+>>>|<row|<cell|<d>+<frac|<y>-<y><rsup|t>|\<gamma\><rsub|t>>-<a><rsub|t>x<rsub|t>-<w>>|<cell|=>|<cell|<0>>>|<row|<cell|<y><rsup|\<top\>><w>>|<cell|=>|<cell|0>>|<row|<cell|x<rsub|t><around*|(|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>-s|)>>|<cell|=>|<cell|0>>|<row|<cell|<around*|(|1-x<rsub|t>|)>s>|<cell|=>|<cell|0>>|<row|<cell|x<rsub|t>,v<rsub|t>,x<rsub|t>,<y>,<w>>|<cell|\<geq\>>|<cell|<0>>>|<row|<cell|x<rsub|t>>|<cell|\<leq\>>|<cell|1,>>>>
  </eqnarray*>

  where the condition with respect to <math|v<rsub|t>> is simplified using
  <math|v<rsub|t>=1-x<rsub|t>> and we have, implicitly, that

  <\equation*>
    <y><rsup|t+1>=<y><rsup|t>-\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)>+\<gamma\><rsub|t><w>
  </equation*>

  for some <math|<w>\<geq\><0>>. Also we notice that due to
  <math|<y><rsup|t+1><rsup|\<top\>><w>=0>, <math|w<rsub|j>=0> whenever
  <math|y<rsub|j><rsup|t+1>\<gtr\>0> and this implies

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>>|<cell|=>|<cell|<around*|\<\|\|\>|<y><rsup|t>-\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)>+\<gamma\><rsub|t><w><rsup|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>-\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>-2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>+\<gamma\><rsub|t><rsup|2><around*|\<\|\|\>|<d>-<a><rsub|t>x<rsub|t>|\<\|\|\>><rsup|2>.>>>>
  </eqnarray*>

  Re-arranging the terms, we obtain

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|2\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>+\<gamma\><rsub|t><rsup|2><around*|\<\|\|\>|<d>-<a><rsub|t>x<rsub|t>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>+\<gamma\><rsub|t><rsup|2>m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>.>>>>
  </eqnarray*>

  Taking <math|\<gamma\><rsub|k>\<equiv\>\<gamma\>>, we take
  <math|<y><rsup|1>=<0>> and telescope from <math|t=1,\<ldots\>,n> to get

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsub|t>|]>>|<cell|\<leq\>>|<cell|<around*|(|2\<gamma\>|)><rsup|-1><big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>><rsup|2>-<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>><rsup|2>+\<gamma\><rsup|2>m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<frac|\<gamma\>|2>m
    n<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>+<around*|(|2\<gamma\>|)><rsup|-1><around*|\<\|\|\>|<y><rsup|1>|\<\|\|\>><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<frac|m
    n<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>|2>\<gamma\>>>>>
  </eqnarray*>

  On the other hand, we know that

  <\equation*>
    <y><rsup|t+1>=<y><rsup|t>-\<gamma\><around*|(|<d>-<a><rsub|t>x<rsub|t>|)>+\<gamma\><w>\<geq\><y><rsup|t>-\<gamma\><around*|(|<d>-<a><rsub|t>x<rsub|t>|)>
  </equation*>

  and\ 

  <\equation*>
    <big|sum><rsub|t=1><rsup|n><a><rsub|t>x<rsub|t>\<leq\><b>+\<gamma\><rsup|-1><big|sum><rsub|t=1><rsup|n><around*|(|<y><rsup|t+1>-<y><rsup|t>|)>\<leq\><b>+\<gamma\><rsup|-1><y><rsup|t+1>,
  </equation*>

  \;

  <\lemma>
    <\em>
      (Informal)

      Let <math|<around*|{|<y><rsup|t>,x<rsub|t>|}>> be a sequence generated
      by the implicit model, then

      <\eqnarray*>
        <tformat|<table|<row|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsub|t>|]>>|<cell|\<leq\>>|<cell|<frac|m
        n<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>|2>\<gamma\>>>|<row|<cell|<big|sum><rsub|t=1><rsup|n><a><rsub|t>x<rsub|t>>|<cell|\<leq\>>|<cell|<b>+\<gamma\><rsup|-1><big|sum><rsub|t=1><rsup|n><around*|(|<y><rsup|t+1>-<y><rsup|t>|)>\<leq\><b>+\<gamma\><rsup|-1><y><rsup|t+1>.>>>>
      </eqnarray*>
    </em>
  </lemma>

  <hrule>

  <with|font-series|bold|Implicit online algorithm>

  Based on <strong|Lemma 1>, we are now ready to propose the implicit online
  algorithn.

  <\eqnarray*>
    <tformat|<table|<row|<cell|<y><rsup|t+1>>|<cell|=>|<cell|arg
    min<rsub|<y>\<geq\><0>> <around*|{|<d><rsup|\<top\>><y>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>|]><rsub|+>+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y>-<y><rsup|t>|\<\|\|\>><rsup|2>|}>>>|<row|<cell|x<rsub|t>>|<cell|=>|<cell|x<around*|(|<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>|]><rsub|+>|)>,>>>>
  </eqnarray*>

  where <math|x<around*|(|<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>|]><rsub|+>|)>>
  denotes the multiplier of linear constraint
  <math|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>\<leq\>s> from the LP
  formulation.

  Before we establish a primal-dual relation for implicit update, we first
  prove the following lemma.

  <\lemma>
    <\em>
      (Bounding the implicit step)

      Given the implicit update of stepsize <math|\<gamma\>>, we have
      <math|<around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|t>|\<\|\|\>>\<leq\>\<gamma\><sqrt|m><around*|(|<wide|a|\<bar\>>+<du>|)>>.
    </em>
  </lemma>

  <\proof>
    Since the implicit model is <math|\<gamma\><rsup|-1>>-strongly convex, we
    invoke three-point lemma and\ 

    <\eqnarray*>
      <tformat|<table|<row|<cell|<d><rsup|\<top\>><y><rsup|t+1>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>+<frac|1|2\<gamma\>><around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|t>|\<\|\|\>><rsup|2>>|<cell|\<leq\>>|<cell|<d><rsup|\<top\>><y><rsup|t>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t>|]><rsub|+>-<frac|1|2\<gamma\>><around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|t>|\<\|\|\>><rsup|2>.>>>>
    </eqnarray*>

    Re-arranging the terms,

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<gamma\><rsup|-1><around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|t>|\<\|\|\>><rsup|2>>|<cell|\<leq\>>|<cell|<d><rsup|\<top\>><around*|(|<y><rsup|t+1>-<y><rsup|t>|)>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t>|]><rsub|+>-<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<sqrt|m><around*|(|<wide|a|\<bar\>>+<du>|)><around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|t>|\<\|\|\>>.>>>>
    </eqnarray*>

    Dividing both sides by <math|\<gamma\><rsup|-1><around*|\<\|\|\>|<y><rsup|t+1>-<y><rsup|t>|\<\|\|\>>>
    completes the proof.
  </proof>

  \;

  Now we are ready to bound the primal-dual gap of implicit update. First we
  have

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<d><rsup|\<top\>><y><rsup|\<ast\>>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|\<ast\>>|]><rsub|+>|]>-\<bbb-E\><rsub|t><around*|[|c<rsub|t>x<rsub|t>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<d><rsup|\<top\>><y><rsup|t>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t>|]><rsub|+>|]>-\<bbb-E\><rsub|t><around*|[|c<rsub|t>x<rsub|t>|]>>>>>
  </eqnarray*>

  and recalling the optimality condition of LP,\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|s>|<cell|\<geq\>>|<cell|<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>>>|<row|<cell|x<rsub|t><around*|(|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>-s|)>>|<cell|=>|<cell|0>>|<row|<cell|<around*|(|1-x<rsub|t>|)>s>|<cell|=>|<cell|0,>>>>
  </eqnarray*>

  we have\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|c<rsub|t>x<rsub|t>>|<cell|=>|<cell|<around*|(|<a><rsub|t><rsup|\<top\>><y><rsup|t+1>+s|)>x<rsub|t>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|)>x<rsub|t>+s
    x<rsub|t>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|)>x<rsub|t>+s>>|<row|<cell|>|<cell|=>|<cell|<around*|(|<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|)>x<rsub|t>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>,>>>>
  </eqnarray*>

  where the last inequality is due to minimization sense of the problem.\ 

  Plugging the relation back,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<d><rsup|\<top\>><y><rsup|t>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t>|]><rsub|+>|]>-<around*|(|<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|)>x<rsub|t>-<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<around*|(|<d>-<a><rsub|t>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>|]>+\<bbb-E\><rsub|t><around*|[|<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t>|]><rsub|+>-<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>|]>+<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<a><rsub|t><rsup|\<top\>><around*|(|<y><rsup|t>-<y><rsup|t+1>|)>x<rsub|t>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<around*|(|<d>-<a><rsub|t><rsup|\<top\>>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>|]>+2n<sqrt|m><wide|a|\<bar\>><around*|\<\|\|\>|<y><rsup|t>-<y><rsup|t+1>|\<\|\|\>>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<around*|(|<d>-<a><rsub|t><rsup|\<top\>>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>|]>+2m
    n\<gamma\><wide|a|\<bar\>><around*|(|<wide|a|\<bar\>>+<du>|)>,>>>>
  </eqnarray*>

  where the second inequality uses\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t>|]><rsub|+>-<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>|]><rsub|+>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\||<a><rsub|t><rsup|\<top\>><around*|(|<y><rsup|t>-<y><rsup|t+1>|)>|\|>\<leq\><around*|\<\|\|\>|<a><rsub|t>|\<\|\|\>>\<cdummy\><around*|\<\|\|\>|<y><rsup|t>-<y><rsup|t+1>|\<\|\|\>>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<sqrt|m><wide|a|\<bar\>><around*|\<\|\|\>|<y><rsup|t>-<y><rsup|t+1>|\<\|\|\>>>>>>
  </eqnarray*>

  and

  <\equation*>
    <a><rsub|t><rsup|\<top\>><around*|(|<y><rsup|t>-<y><rsup|t+1>|)>x<rsub|t>\<leq\><around*|\||<a><rsub|t><rsup|\<top\>><around*|(|<y><rsup|t>-<y><rsup|t+1>|)>|\|>\<cdummy\><around*|\||x<rsub|t>|\|>\<leq\><sqrt|m><wide|a|\<bar\>><around*|\<\|\|\>|<y><rsup|t>-<y><rsup|t+1>|\<\|\|\>>.
  </equation*>

  The last inequality invokes <strong|Lemma 2> bounding
  <math|<around*|\<\|\|\>|<y><rsup|t>-<y><rsup|t+1>|\<\|\|\>>\<leq\>\<gamma\><sqrt|m><around*|(|<wide|a|\<bar\>>+<du>|)>>.

  <hrule>

  <with|font-series|bold|Stability of implicit model>

  Before arriving at our final results, we bound the iteration
  <math|<around*|\<\|\|\>|<y><rsup|t>|\<\|\|\>>> using a stability argument
  of stochastic proximal point method [1],\ 

  <\lemma>
    <\em>
      (Stability of stochastic proximal point)

      Let <math|<around*|{|<y><rsup|t>|}>> be the sequence of iterates
      generated by implicit model, then

      <\equation*>
        \<bbb-E\><rsub|t><around*|[|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>>|]>\<leq\><sqrt|Y<rsub|s>+m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsup|2>t>+Y<rsub|l>\<leq\><sqrt|<around*|(|<frac|<wide|c|\<bar\>>|<dl>>|)><rsup|2>+
        m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsup|2>t>+<frac|<wide|c|\<bar\>>|<dl>>,
      </equation*>

      where <math|Y<rsub|s>\<assign\>inf<rsub|<y>\<in\>\<cal-Y\><rsup|*\<ast\>>><around*|\<\|\|\>|<y>|\<\|\|\>>>
      and<space|1em><math|Y<rsub|l>\<assign\>sup<rsub|<y>\<in\>\<cal-Y\><rsup|*\<ast\>>><around*|\<\|\|\>|<y>|\<\|\|\>>>
      denote the minimum/maximum norm of dual optimal solutions.
    </em>
  </lemma>

  <\proof>
    Let <math|\<cal-Y\><rsup|\<ast\>>> denote the set of dual optimal
    solutions. Then we invoke <strong|Corollary 3.1> from [1] and
    successively deduce that

    <\eqnarray*>
      <tformat|<table|<row|<cell|>|<cell|>|<cell|\<bbb-E\><rsub|t><around*|[|<around*|\<\|\|\>|<y><rsup|t+1>-Proj<rsub|\<cal-Y\><rsup|\<ast\>>><around*|(|<y><rsup|t+1>|)>|\<\|\|\>><rsup|2>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|\<\|\|\>|<y><rsup|1>-Proj<rsub|\<cal-Y\><rsup|\<ast\>>><around*|(|<y><rsup|1>|)>|\<\|\|\>><rsup|2>+m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsup|2>t>>|<row|<cell|>|<cell|=>|<cell|<around*|\<\|\|\>|<y><rsup|\<ast\>>|\<\|\|\>><rsup|2>+m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsup|2>t.>>>>
    </eqnarray*>

    By Jensen's inequality,

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<bbb-E\><rsub|t><around*|[|<around*|\<\|\|\>|<y><rsup|t+1>|\<\|\|\>>|]>>|<cell|\<leq\>>|<cell|\<bbb-E\><rsub|t><around*|[|<around*|\<\|\|\>|<y><rsup|t+1>-Proj<rsub|\<cal-Y\><rsup|\<ast\>>><around*|(|<y><rsup|t+1>|)>|\<\|\|\>>|]>+\<bbb-E\><rsub|t><around*|[|<around*|\<\|\|\>|Proj<rsub|\<cal-Y\><rsup|\<ast\>>><around*|(|<y><rsup|t+1>|)>|\<\|\|\>>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<sqrt|<around*|\<\|\|\>|<y><rsup|\<ast\>>|\<\|\|\>><rsup|2>+m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsup|2>t>+<frac|<wide|c|\<bar\>>|<dl>>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<sqrt|<around*|(|<frac|<wide|c|\<bar\>>|<dl>>|)><rsup|2>+m<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>\<gamma\><rsup|2>t>+<frac|<wide|c|\<bar\>>|<dl>>>>>>
    </eqnarray*>

    where the first relation uses triangle inequality and the second uses
    <math|<around*|\<\|\|\>|<y><rsup|\<ast\>>|\<\|\|\>>\<leq\><around*|\<\|\|\>|<y><rsup|\<ast\>>|\<\|\|\>><rsub|1>\<leq\><frac|<wide|c|\<bar\>>|<dl>>>.
  </proof>

  <\remark>
    We note that compared to the subgradient case, <math|Y<rsub|s>,Y<rsub|l>>
    may be significantly smaller than <math|<frac|<wide|c|\<bar\>>|<dl>>>,
    thus enhancing the practical performance of the implicit model.
  </remark>

  [1] Asi, H. & Duchi, J. C. Stochastic (Approximate) Proximal Point Methods:
  Convergence, Optimality, and Adaptivity.<nbsp><with|font-shape|italic|Arxiv><nbsp>(2018)
  doi:10.48550/arxiv.1810.05633.

  <hrule>

  <with|font-series|bold|Balancing regret and violation>.

  As we do for subgradient method, we now have

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><around*|[|R<rsub|n><rsup|\<ast\>>-R<rsub|n>|]>>|<cell|\<leq\>>|<cell|<big|sum><rsub|t=1><rsup|n>\<bbb-E\><rsub|t><around*|[|<around*|(|<d>-<a><rsub|t><rsup|\<top\>>x<rsub|t>|)><rsup|\<top\>><y><rsup|t>|]>+2m
    n\<gamma\><wide|a|\<bar\>><around*|(|<wide|a|\<bar\>>+<du>|)>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<around*|[|<frac|<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>|2>+2<wide|a|\<bar\>><around*|(|<wide|a|\<bar\>>+<du>|)>|]>m
    n\<gamma\>>>>>
  </eqnarray*>

  and\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><around*|[|<around*|\<\|\|\>|<A><x>-<b>|\<\|\|\>><rsub|+>|]>>|<cell|\<leq\>>|<cell|\<gamma\><rsup|-1>\<bbb-E\><rsub|t><around*|[|<around*|\<\|\|\>|<y><rsup|n+1>|\<\|\|\>>|]>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<sqrt|<around*|(|<frac|<wide|c|\<bar\>>|<dl>>|)><rsup|2>\<gamma\><rsup|-2>+m
    n<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2>>+<frac|<wide|c|\<bar\>>|\<gamma\><dl>>>>|<row|<cell|>|<cell|\<leq\>>|<cell|<sqrt|m
    n><around*|(|<wide|a|\<bar\>>+<du>|)>+<frac|2<wide|c|\<bar\>>|\<gamma\><dl>>,>>>>
  </eqnarray*>

  where we use the relation <math|<sqrt|x+y>\<leq\><sqrt|x>+<sqrt|y>>.

  Taking a trade-off between the two terms, we have\ 

  <\equation*>
    \<gamma\><rsup|\<ast\>>=<frac|1|<sqrt|m
    n>><around*|[|<frac|<around*|(|<wide|a|\<bar\>>+<du>|)><rsup|2><dl>|2<wide|c|\<bar\>>>+<frac|2<dl>|<wide|c|\<bar\>>><wide|a|\<bar\>><around*|(|<wide|a|\<bar\>>+<du>|)>|]><rsup|-1/2>
  </equation*>

  and the optimal constraint violation and regret are both
  <math|\<cal-O\><around*|(|<sqrt|m n>|)>>.

  <hrule>

  <with|font-series|bold|Solving the proximal point subproblem>

  Last we consider the practical solution of the proximal point subproblem.

  <\eqnarray*>
    <tformat|<table|<row|<cell|<y><rsup|t+1>>|<cell|=>|<cell|arg
    min<rsub|<y>\<geq\><0>> <around*|{|<d><rsup|\<top\>><y>+<around*|[|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>|]><rsub|+>+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y>-<y><rsup|t>|\<\|\|\>><rsup|2>|}>>>>>
  </eqnarray*>

  We do case analysis on <math|<y><rsup|t+1>>.

  <with|font-series|bold|Case 1.> <math|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>\<gtr\>0>

  In this case we have <math|<y><rsup|t+1>=<y><rsup|t>-\<gamma\><rsub|t><around*|(|<d>-<a><rsub|t>|)>>

  <with|font-series|bold|Case 2.> <math|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>\<less\>0>

  In this case we have <math|<y><rsup|t+1>=<y><rsup|t>-\<gamma\><rsub|t><d>>

  <with|font-series|bold|Case 3.> <math|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y><rsup|t+1>=0>

  This case implies an optimization problem

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|<y>>>|<cell|<d><rsup|\<top\>><y>+<frac|1|2\<gamma\><rsub|t>><around*|\<\|\|\>|<y>-<y><rsup|t>|\<\|\|\>><rsup|2>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|c<rsub|t>-<a><rsub|t><rsup|\<top\>><y>=0>|<cell|>>|<row|<cell|>|<cell|<y>\<geq\><0>.>|<cell|>>>>
  </eqnarray*>
</body>

<\initial>
  <\collection>
    <associate|font-base-size|10>
    <associate|page-medium|paper>
  </collection>
</initial>