#lipidomics analysis using lipidr 
#the data used is the sample data from the lipidr website 


#run as administrator 
library(devtools)   
install_github("ahmohamed/lipidr")

dm_path = read.csv("C:\\Users\\Aashi\\Desktop\\lipidomics_analyis\\test_data_matrix.csv")
ta_path = read.csv("C:\\Users\\Aashi\\Desktop\\lipidomics_analyis\\test_annotation_data.csv")
d <- as_lipidomics_experiment(read.csv("C:\\Users\\Aashi\\Desktop\\lipidomics_analyis\\test_data_matrix.csv"))
d <- add_sample_annotation(d, ta_path)
# visualize lipid intensity across the samples 
# https://www.lipidr.org/articles/examples/mw_integration.html
plot_samples(d, type="tic", log=TRUE)

d <- set_logged(d, "Area", TRUE)
d <- set_normalized(d, "Area", TRUE)
plot_samples(d, "boxplot")
mvaresults = mva(d, measure="Area", method="PCA")
plot_mva(mvaresults, color_by="SampleType", components = c(1,2))

two_group <- de_analysis(d, Cancer-Benign, Cancer-Metastasis)
plot_results_volcano(two_group)