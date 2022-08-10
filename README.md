# RingLWE
* Chi2Attack.sage

This is the source code for chi-square attack proposed at「Attacks on the Search-RLWE problem with small errors」.
Note that even small parameters can be very time consuming.


* Prime_residue_degree_Chi2Attack_Basic.sage, Prime_residue_degree_Chi2Attack_Improved.sage
attack using cosets $\mathbb{F}_{q^f}/\mathbb{F}_q$.


* 4_residue_degree_Chi2Attack_Basic.sage, 4_residue_degree_Chi2Attack_Improved.sage
This is the source code especially for f=4. In this attack, we using cosets $\mathbb{F}_{q^4}/\mathbb{F}_{q^2}$.

* 6_residue_degree_Chi2Attack_Basic.sage, 6_residue_degree_Chi2Attack_Improved.sage
This is the source code especially for f=6. In this attack, we using cosets $\mathbb{F}_{q^6}/\mathbb{F}_{q^3}$.

* DirectCycSampler.sage
* ExtendCyclotomic.sage
* mega.sage
* misc.sage
* MyLatticeSamplar.sage
* PeikertSampler.sage
* QuadCyclotomic.sage
* subcycsampler.sage
* SubgroupModm.sage
in the "/sampling" are the code published in HaoChen's「Attacks on the Search-RLWE problem with small errors」.
