<TeXmacs|2.1.1>

<style|generic>

<\body>
  <doc-data|<doc-title|Fast initialization for linear programs>>

  In this note we discuss fast algorithm to start LPs

  <\eqnarray*>
    <tformat|<table|<row|<cell|max<rsub|x>>|<cell|c<rsup|\<top\>>x>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|A x=b>|<cell|>>|<row|<cell|>|<cell|0\<leq\>x\<leq\>e>|<cell|>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|y>>|<cell|b<rsup|\<top\>>y+e<rsup|\<top\>>s>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|s\<geq\>c-A<rsup|\<top\>>y>|<cell|>>|<row|<cell|>|<cell|y,s\<geq\>0>|<cell|>>>>
  </eqnarray*>

  Using reformulation

  <\eqnarray*>
    <tformat|<table|<row|<cell|min<rsub|y>>|<cell|<around*|(|<frac|b|n>|)><rsup|\<top\>>y+<frac|1|n><big|sum><rsub|i=1><rsup|n><around*|[|c<rsub|i>-a<rsub|i><rsup|\<top\>>y|]><rsub|+>>|<cell|>>|<row|<cell|<text|subject
    to>>|<cell|y\<geq\>0.>|<cell|>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|d<rsup|\<top\>>y+<frac|1|n><big|sum><rsub|i=1><rsup|n><around*|[|c<rsub|i>-a<rsub|i><rsup|\<top\>>y|]><rsub|+>-d<rsup|\<top\>>y<rsup|\<ast\>>+<frac|1|n><big|sum><rsub|i=1><rsup|n><around*|[|c<rsub|i>-a<rsub|i><rsup|\<top\>>y<rsup|\<ast\>>|]><rsub|+>>>|<row|<cell|>|<cell|\<geq\>>|<cell|<frac|1|n><big|sum><rsub|i=1><rsup|n>max<around*|{|<around*|(|a<rsub|i>+b|)><rsup|\<top\>><around*|(|y-y<rsup|\<ast\>>|)>,b<rsup|\<top\>><around*|(|y-y<rsup|\<ast\>>|)>|}>>>|<row|<cell|>|<cell|\<geq\>>|<cell|<around*|\<\|\|\>|y-y<rsup|\<ast\>>|\<\|\|\>>>>>>
  </eqnarray*>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
  </collection>
</initial>