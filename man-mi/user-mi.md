User Guide for ChaoMI
================
<font size="4">
* [Overview](#view2)
  * [How to Cite](#cite2)
* [Data Settings](#data-settings2)
  * [Dataset](#dataset2)
  * [Demo Data](#demo2)
  * [Upload Data](#upload2)
* [General Settings](#gen-settings2)
  * [Method](#method2)
  * [Bootstraps](#bootstraps2)
  * [Confidence level](#confidence2)
* [Output](#output2)
  * [Data Summary](#summary2)
  * [Esitmation](#chaomi2)
  * [Visualization](#vis2)
* [References](#reference2)

* * * * * * * *

<h2 id="view2">Overview</h2>
The program <strong>`ChaoMI`</strong> ( <strong>`Mutual Information`</strong> proposed by <strong>`Chao`</strong> et al.) is written in the <a href="http://www.r-project.org/" target="_blank">R</a> language and the interactive web application is built by using <a href="http://www.rstudio.com/shiny" target="_blank">Shiny</a>. The user provides a 
matrix of two-dimension data, in which cells with positive integers indicate the frequency of interaction between a pair of species (variables). <strong>`ChaoMI`</strong> computes the Mutual Information estimators, bootstrap standard error and confidence intervals.

_<h3 id="cite2">How to Cite</h3>_
<font color="ff0000">If you use <strong>`ChaoEntropy Online`</strong> to obtain results for publication, you should cite the papers (Chao, A., Wang, Y. T. and Jost, L. (2013) Entropy and the species accumulation curve: a novel entropy estimator via discovery rates of new species. _Methods in Ecology and Evolution_, __4__, 1091-1100.) along with the following reference for `ChaoEntropy Online`:

<p style="padding-left: 30px;"> Chao, A., Lee, Y.-H. and Tseng, K.-S. (2014). ChaoEntropy Online.</font>
  
To help refine <strong>`ChaoEntropy Online`</strong>, your comments or feedbacks would be welcome (please send them to chao@stat.nthu.edu.tw).

<h2 id="data-settings2">Data Settings</h2>

_<h3 id="dataset2">Dataset</h3>_
Some demonstration datasets are used for illustration. 

* Schleuning et al. (2010): This describes a plant-frugivore network from a Kenyan rain forest. Interaction frequency in the quantitative interaction matrix was defined as the number of fruit-eating individuals on a plant species independent of fruit handling (including swallowing, pecking, and dropping fruits).
  For more detail, you can see the source of this data : http://www.nceas.ucsb.edu/interactionweb/html/schleuning-et-al-2010.html
  
* Bluthgen et al. (2004): This describes de ant-plant interaction network in the rainforest at the Australian Canopy Crane in Cape Tribulation, Far North Queensland, Australia (16°07' S, 145°27' E, 80 m a.s.l.). Frequency of interaction was estimated as the number of ant individuals of a particular species collected on a plant species. Data are presented as an interaction frequency matrix, in which cells with positive integers indicate the frequency of interaction between a pair of species, and cells with zeros indicate no interaction.
  For more detail, you can see the source of this data : http://www.nceas.ucsb.edu/interactionweb/html/bluthgen_2004.html

_<h3 id="demo2">Demo Data</h3>_
<strong>`ChaoMI`</strong> provides three demo dataset from Schleuning et al. (2010).
We according to forest strata deviding the data into three class:
  1. `canopy`
  2. `midstory`
  3. `understory`

You can compare its interaction network from different forest strata.

_<h3 id="upload2">Upload Data</h3>_
<strong>`ChaoMI`</strong> also provides a upload function. Click the **Upload data** button, there will show the **file choose** button to upload the data. **User's data MUST be .csv**. After upload the data, the file name which your uploaded will be listed in the box. 

We provides three dataset to download from the Bluthgen et al. (2004), according to the subfamily of ants deviding the data into three group:
  1. <a href="https://dl.dropboxusercontent.com/s/jgodzi1orkfk746/Dolichoderinae.csv">`Dolichoderinae`(Click to download)</a>
  2. <a href="https://dl.dropboxusercontent.com/s/77c2pst6z9l4bq5/Formicinae.csv">`Formicinae`(Click to download)</a>
  3. <a href="https://dl.dropboxusercontent.com/s/dtn73daay1j6jp1/Myrmicinae.csv">`Myrmicinae`(Click to download)</a>

<font color="ff0000">Notice</font>: In the upload files, the variable names are better as simple as possible. It can't include symbol like <font color="ff0000">`' '`</font>, <font color="ff0000">`" "`</font>. If the Data Viewer was blank after click Run button, check the variable names and see the demo download files above.

<h2 id="gen-settings2">General settings</h2>
_<h3 id="method2">Method</h3>_
This is a checkbox for select the method which constructed to estimate mutual information.
  * `Chao` estimator, for detail see reference 1.
  * `ChaoShen` estimator, for detail see reference 3.
  * `Jackknife` estimator, for detail see reference 5.
  * `Bias Correct 1` estimator.
  * `Bias Correct 2` estimator.
  * `Observed` estimator.

_<h3 id="bootstraps2">Bootstraps</h3>_
Number of bootstraps (say B) is an integer specifying the number of replications for bootstrap resampling scheme in computing variance. Refer to Chao et al. (2013) for details. Default is `100`. To save running time, we recommend that 100 or 200 bootstraps will be sufficient for most applications.  

_<h3 id="confidence2">Confidence level</h3>_
The confidence level is a positive number is less than or equal to 1. The default is `0.95`.

<h2 id="output2">Output</h2>
_<h3 id="summary2">Data Summary</h3>_
This tab panel shows basic data information for the selected data. The output variables are interpreted at the first column.

_<h3 id="chaomi2">Estimation</h3>_
This tab panel shows the main output for <strong>`ChaoMI`</strong>. It show a table of various mutual information estimators, their standard error and confidence interval which the method you choose. You also can click [Download as csv file]() to download the output table.

_<h3 id="vis2">Visualization</h3>_
This tab panel shows the interactive estimator and confidence interval plot.

<h2 id="reference2">References</h2>
1. Chao, A., Wang, Y. T. and Jost, L. (2013) Entropy and the species accumulation curve: a novel entropy estimator via discovery rates of new species. _Methods in Ecology and Evolution_, __4__, 1091-1100.

2. Chao, A. & Jost, L. (2012) Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. _Ecology_, __93__, 2533-2547.

3. Chao, A. & Shen, T.J. (2003) Nonparametric estimation of Shannon's index of diversity when there are unseen species. _Environmental and Ecological Statistics_, __10__, 429-443.

4. Zahl, S. (1977) Jackknifing an index of diversity. _Ecology_, __58__, 907-913.

5. Blüthgen, N., N. E. Stork, and K. Fiedler. 2004. Bottom-up control and co-occurrence in complex communities: honeydew and nectar determine a rainforest ant mosaic. _Oikos_, __106__, 344-358.

6. Matthias Schleuning, Nico Blüthgen, Martina Flörchinger, Julius Braun, H. Martin Schaefer, Katrin Böhning-Gaese, (2011). Specialization and interaction strength in a tropical plant-frugivore network differ among forest strata. _Ecology_, __92(1)__, 26-36.
</font>