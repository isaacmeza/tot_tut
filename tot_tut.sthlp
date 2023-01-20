{smcl}
{* *!version 1.0.3  2023-01-19}{...}
{viewerjumpto "Syntax" "tot_tut##syntax"}{...}
{viewerjumpto "Description" "tot_tut##description"}{...}
{viewerjumpto "Options" "tot_tut##options"}{...}
{viewerjumpto "Examples" "tot_tut##examples"}{...}
{viewerjumpto "Stored results" "tot_tut##stored_results"}{...}
{viewerjumpto "References" "tot_tut##references"}{...}
{viewerjumpto "Authors" "tot_tut##authors"}{...}

{title:Title}

{p 4 8}{cmd:tot_tut} {hline 2} Implementation for the estimation of treatment on the treated (ToT), treatment on the untreated (TuT) and the average treatment effect (ATE) jointly using the design introduced in "The limits of self-commitment and private paternalism". Other causal statistics of interest like the average selection bias and the average selection on the level are also provided. {p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:tot_tut} {it:depvar} {it:treatvar} {it:choicevar} {ifin} 
[{cmd:,} 
{cmd:vce(robust | cluster {it:clustvar}{cmd:)}} 
]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:tot_tut} estimates jointly the treatment on the treated, treatment on the untreated, the average treatment effect, selection on gains, selection bias, and selection on the level, exploiting a design with three arms: a control arm, a forced arm and a choice arm. The specification strategy involves estimating two iv regressions per each selection estimand, and jointly obtaining errors. Details on the implementation can be found  in Section 5 and Appendix of the paper  {browse "https://isaacmeza.github.io/personal//files/donde.pdf": "The limits of self-commitment and private paternalism"}. 

Note : Jointly inference for selection on gains, selection bias, and selection on the level is not provided. 


{marker arguments}{...}
{title:Arguments}
{dlgtab:Arguments}

{p 8 8} {it:depvar}, this is the outcome of interest.{p_end}

{p 8 8} {it:treatvar}, categorical variable indicating treatment status: control arm (0), forced arm (1), choice arm (2). {p_end}

{p 8 8} {it:choicevar}, binary variable indicating choice.{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Options}

{p 4 8}{cmd:vce(robust | cluster {it:clustvar}{cmd:)}} specifies the type of standard error reported, which includes types that are robust to some kinds of misspecification (robust - the default), and that allow for intragroup correlation (cluster {it:clustvar}).{p_end}

{hline}


{marker examples}{...}
{title: Examples}

{title: "The limits of self-commitment and private paternalism"}

{p 4 8}Setup{p_end}
{p 8 8}{stata "use tot_tut_commitment.dta, clear"}{p_end}
{p 8 8}{stata "gen x0 = -(Z==2)*(choose==0)"}{p_end}
{p 8 8}{stata "gen x1 = (Z==2)*(choose==1)"}{p_end}
{p 8 8}{stata "gen z0_ = -(Z==0)"}{p_end}
{p 8 8}{stata "gen z0 = (Z==0)"}{p_end}
{p 8 8}{stata "gen z1 = (Z==1)"}{p_end}
{p 8 8}{stata "gen z2 = (Z==2)"}{p_end}

{p 4 8}ToT & ATE using ivregress {p_end}
{p 8 8}{stata "ivregress 2sls apr z1 (x1 = z2), vce(cluster clustvar)"}{p_end}

{p 4 8}TuT & ATE using ivregress {p_end}
{p 8 8}{stata "ivregress 2sls apr z0_ (x0 = z2), vce(cluster clustvar)"}{p_end}

{p 4 8}Simultaneous inference for ToT & TuT - selection on gains{p_end}
{p 8 8}{stata "tot_tut apr Z choose, vce(cluster clustvar)"}{p_end}


{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:tot_tut} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}

{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations.{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom.{p_end}
  
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient fector.{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators.{p_end}

{marker references}{...}
{title:References}

{p 4 8}{browse "https://isaacmeza.github.io/personal//files/donde.pdf": DiTraglia, McIntosh, Meza, Seira, Sadka.} "The limits of self-commitment and private paternalism". Working paper. {p_end}


{marker authors}{...}
{title:Authors}

{p 4 8}Meza Lopez Isaac; Harvard, Economics.
{browse "mailto:isaacmezalopez@g.harvard.edu":isaacmezalopez@g.harvard.edu}.{p_end}