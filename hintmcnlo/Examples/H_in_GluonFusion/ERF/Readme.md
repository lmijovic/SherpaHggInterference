# Changes wrt YR4  

Need to adapt selectors and isolation: 

```
AddOns/Higgs/Higgs_Selector.C 
str<<"HiggsFinder pt1 pt2 eta mmin mmax [dr epspt]";

```

```
IsolationCut 22 &lt;dR&gt; &lt;exponent&gt; &lt;epsilon&gt;
implements the Frixione isolation cone [Frixione1998]

```

```
(selector){
  HiggsFinder 5 5 2.5 100 150;
  IsolationCut 22 0.1 2 0.025;
}(selector); 

```

