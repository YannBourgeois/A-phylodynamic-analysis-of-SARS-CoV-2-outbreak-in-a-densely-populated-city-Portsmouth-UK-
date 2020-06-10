# A phylodynamic analysis of SARS-CoV-2 outbreak in a densely populated city (Portsmouth, UK).
A repository with reports on phylodynamic analyses realized on COVID-19 sequences collected at the University of Portsmouth/Portsmouth Hospital

Yann Bourgeois (1), Chauhan Anoop (2), Angela Beckett (1), Glaysher Sharon (2), Elliott Scott (2), Ethan Butcher (2), Sarah Wyllie (2), Lloyd Allyson (2), Garry Scarlett (1), Bicknell Kelly (2), Samuel Robson (1,3)

1.	University of Portsmouth
2.	Portsmouth Hospitals NHS Trust
3.	Centre for Enzyme Innovation, University of Portsmouth

**Summary**

Tracking epidemics requires heavy logistics, and may miss many cases due to insufficient testing capabilities or fast spread. The examination of viral genomes may help from this perspective. Genealogies of viral sequences contain information about the processes that generate them. Recent advances in the modelling of the coalescence process during epidemics give access to relevant epidemiological parameters such as the effective reproductive number or prevalence. In this report, we analysed SARS-CoV-2 sequences collected at Portsmouth Hospitals. From an epidemiological perspective, Portsmouth is an important city in the UK, being second only to London for its population density. It is therefore vulnerable to a fast spread of coronavirus. We use the birth-death coalescent framework implemented in two BEAST2 packages to reconstruct the recent dynamic of the epidemic and provide a proof of concept for further, more detailed analyses. 

**Methods**

_Data filtering_
We used a set of 105 genomes sequenced using a MinION procedure (see the ARTIC network protocol (https://artic.network/ncov-2019)). We aligned consensus sequences on the reference genome (GENBANK accession MN908947.3) using MAFFTv7.453 (Katoh and Standley, 2013). We masked the alignment to exclude spurious variants due to poor mapping and sites prone to homoplasy (http://virological.org/t/issues-with-sars-cov-2-sequencing-data/473). Masking was performed using the maskfasta option in BEDTOOLS (Quinlan and Hall, 2010). After masking, we excluded sequences with more than 15% missing data, estimating the number of ambiguous sites with the nuc option in BEDTOOLS. The birth-death models underlying these two analyses assume a uniform sampling from the epidemic. This means that oversampling of a local transmission cluster may bias the results. To limit overrepresentation of local transmission cases, we removed identical sequences from the dataset, using a UPGMA distance matrix estimated in MEGA X (Kumar et al., 2018). We only kept sequences from patients living in Portsmouth region. The final dataset contained 35 sequences sampled between April 27th and May 25th.

_Birth-death coalescent framework._
We used two different Bayesian methods, implemented in the BEAST2 software (Drummond and Rambaut, 2007). The first package, BDSKY (Stadler et al., 2013), applies the coalescent skyline plot framework to viral sequences. It can estimate the sampling rate, the rate at which individuals become non-infectious, or the effective reproduction number (_R_). All these parameters can change through time in a piecewise way. The second package, EpiInf (Vaughan et al., 2019), is based on the same framework, but also attempts are estimating the evolution of prevalence through time. We used a SIR model, where individuals that survive infection (“I”) are removed (“R”) from the susceptive pool (“S”) and cannot be infected a second time. We adopted this analytical protocol given its ability to reconstruct the dynamics of COVID-19 infections for several outbreaks (http://virological.org/t/phylodynamic-analyses-of-outbreaks-in-china-italy-washington-state-usa-and-the-diamond-princess/439).

_Model of sequence evolution_
We used a HKY model of sequence evolution (Hasegawa et al., 1985), and fixed the substitution rate to 0.001 substitutions/site/year (Duchene et al., 2020). Note that parameters such as prevalence or times will vary according to this rate: using a lower rate would increase estimates. Given that other commonly used estimates of the substitution rate are lower (see for example http://virological.org/t/update-2-evolutionary-epidemiological-analysis-of-128-genomes/423), it is likely that our inference does not overestimate times or prevalence.

_Priors_
For the first analysis (BDSKY), we used a multi-rho sampling model, which takes into account the fact that we sampled multiple sequences simultaneously. We allowed the effective reproduction number (R) to change four times after the start of the epidemics. Given the relatively low number of sequences at this stage, we keep all other parameters constant over time. Priors for parameters are given in Table 2. We used a fixed removal rate of 36.5 year-1 (365/10). Removal rate is the inverse of the time (in years) for an individual to become non-infectious, which is considered reasonable (see http://virological.org/t/update-2-evolutionary-epidemiological-analysis-of-128-genomes/423, (Li et al., 2020)).
For the second analysis (EpiInf), we specified that samples taken from patients at the hospital coincided with a full removal from the infectious pool (Removal proportion of 1.0). Due to poor mixing of the chains, we did not estimate the initial susceptible population size, and set it to 100,000, which is in the range of the local demographic count.


Model | Parameter | Prior
-- | -- | --
Core | Substitution rate (strict   clock) | Fixed at 0.001   sub/site/year
-- | Evolution model | HKY. Kappa ~ lognormal(1.0,   1.25)
BDSKY/EpiInf | Reproductive numbers (x4 in   BDSKY) | ~LogNormal(0,1)
-- |Recovery Rate (=Becoming   Uninfectious) | Fixed at 36.5
-- |Sampling proportion | ~ 1/X, within interval (0,   1)
-- |Infection rate (=birth   rate) | ~ 1/X
-- |Initial susceptible   population size | Fixed at 100,000
-- |Epidemic origin time | ~ 1/X, within interval (0,   20)

_Table 1. Priors used for BDSKY and EpiInf analyses in BEAST2_

_Markov Chain Monte Carlo (MCMC) procedure and replication_
We ran ten independent chains for both analyses, each lasting two million generations. States were sampled every 1000 generations. For EpiInf, which relies on a particle filter algorithm for estimating prevalence trajectories (see (Vaughan et al., 2019)), we used the default 100 particles per tree prior evaluation. Runs were combined together to obtain posterior distributions for all parameters, removing the first 500,000 generations as burn-in. Visual examination of posterior traces was performed using Tracer v1.7.1 to ensure convergence and good mixing of chains (Rambaut et al., 2018). 

**Results**

_Origin of the epidemic_
Using BDSKY, we could estimate the time at which the epidemic most likely started. We obtained a median estimate of 0.32 years before last sampling, with a 95% highest posterior density (HPD) of [0.27, 0.37]. This assumes a substitution rate of 0.001 substitution/year, and places the start of the epidemic in Portsmouth’s region between January 11th and February 17th.
Using EpiInf, we obtained a median estimate of 0.21 years, with a broader highest posterior density of [0.16, 0.27], placing the start of the epidemic between February 17th and March 28th. 
Evolution of _R_, prevalence and slowing down of the epidemic
The four estimates of _R_ since the start of the epidemic are shown on Figure 1. There is a slow increase of _R_ during the first months of the epidemic, with a steep increase in April, followed by a strong decrease suggesting a negative growth rate of the virus in the first weeks of May. Again, we remind that the exact times may vary due to uncertainties on the substitution rate and sampling.

![image](https://user-images.githubusercontent.com/6105724/84243321-eee84700-aaf9-11ea-816a-d0acb88a9a61.png)

_Figure 1. Estimates of the effective reproductive number R over five distinct time intervals since the start of the epidemic in Portsmouth. Values above 1 indicate growth in the number of cases, while values below 1 suggest a decline._

Results from BDSKY were nonetheless qualitatively consistent with estimates of prevalence inferred by EpiInf (Figure 2). There was a clear rise in prevalence in the last weeks of April, followed by a slow-down. Median estimate of prevalence at the time of sampling (May 25th) was 93, with a HPD of [61, 129].

![image](https://user-images.githubusercontent.com/6105724/84243526-3cfd4a80-aafa-11ea-8759-92fe35e802d3.png)

_Figure 2 Estimates of prevalence (number of cases) through time obtained from EpiInf analyses. Red lines correspond to different samples from the posterior set of prevalence trajectories._

**Discussion and future directions**

This preliminary analysis will require further sampling across Portsmouth and Hampshire to ensure the robustness of results and improve our estimates of prevalence and origin of the outbreak. The dynamic that we highlight is nevertheless consistent with observations in the field, which reported a steep increase in the number of cases at Portsmouth Hospital in the last weeks of April. The high population density in Portsmouth did not seem to lead to substantial increases in the reproductive number at this stage, which remained at its maximum in the range of other studies (R between 3 and 5).
Quantifying the impact of COVID-19 is of primordial importance to design public policies that protect the most vulnerable populations. Metadata on patients’ geographic origin, socio-economic background and age will be primordial from this perspective, and are currently collected. We will use other methods of phylodynamic inference that can handle multiple reproductive numbers and propagation rates in connected subpopulations (BDMM, (Kühnert et al., 2016)). These results may help to estimate how the virus spread in different subpopulations, and inform public policies in case of a second wave of infection.

**References**

Drummond AJ, Rambaut A (2007). BEAST : Bayesian evolutionary analysis by sampling trees. 8: 1–8.

Duchene S, Featherstone L, Haritopoulou-sinanidou M, Rambaut A, Baele G (2020). Temporal signal and the phylodynamic threshold of SARS-CoV-2. bioRxiv.

Hasegawa M, Kishino H, Yano T aki (1985). Dating of the human-ape splitting by a molecular clock of mitochondrial DNA. J Mol Evol 22: 160–174.

Katoh K, Standley DM (2013). MAFFT multiple sequence alignment software version 7: Improvements in performance and usability. Mol Biol Evol 30: 772–780.

Kühnert D, Stadler T, Vaughan TG, Drummond AJ (2016). Phylodynamics with Migration: A Computational Framework to Quantify Population Structure from Genomic Data. Mol Biol Evol 33: 2102–2116.

Kumar S, Stecher G, Li M, Knyaz C, Tamura K (2018). MEGA X: Molecular evolutionary genetics analysis across computing platforms. Mol Biol Evol 35: 1547–1549.

Li Q, Guan X, Wu P, Wang X, Zhou L, Tong Y, et al. (2020). Early transmission dynamics in Wuhan, China, of novel coronavirus-infected pneumonia. N Engl J Med 382: 1199–1207.

Quinlan AR, Hall IM (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26: 841–842.
Rambaut A, Drummond AJ, Xie D, Baele G, Suchard MA (2018). Posterior summarization in Bayesian phylogenetics using Tracer 1.7. Syst Biol 67: 901–904.

Stadler T, Kühnert D, Bonhoeffer S, Drummond AJ (2013). Birth-death skyline plot reveals temporal changes of epidemic spread in HIV and hepatitis C virus (HCV). Proc Natl Acad Sci U S A 110: 228–233.

Vaughan TG, Leventhal GE, Rasmussen DA, Drummond AJ, Welch D, Stadler T, et al. (2019). Estimating Epidemic Incidence and Prevalence from Genomic Data. Mol Biol Evol 36: 1804–1816.

