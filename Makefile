num_cores = 1

experiments = ber_id ber_eq ber_toe poi_id bin_id

all: $(experiments) sims_plots real_data

clean:
	rm -r output

$(experiments): %:  output/%/setup.rda output/%/gt.qs output/%/mle.qs output/%/bt.qs output/%/armb.qs output/%/blb.qs

output/%/setup.rda: R/sett/%.R R/sett/common.R R/sett/common_ri.R
	@echo "\n*********************************\
				 \n* $*: Initialisation\
				 \n*********************************\n"
	mkdir -p output/$*/
	Rscript R/sett/common.R
	Rscript R/sett/common_ri.R
	Rscript R/sett/$*.R

output/%/gt.qs: R/ground_truth.R output/%/setup.rda
	@echo "\n****************************************************\
				 \n* $*: Computing the ground truth on $(num_cores) cores \
				 \n****************************************************\n"
	Rscript R/ground_truth.R $* $(num_cores)

output/%/mle.qs: R/mle.R output/%/setup.rda
	@echo "\n**********************************************\
				 \n* $*: Computing the mle on $(num_cores) cores \
				 \n**********************************************\n"
	Rscript R/mle.R $* $(num_cores)

output/%/bt.qs: R/bootstrap.R output/%/setup.rda output/%/mle.qs
	@echo "\n**********************************************\
				 \n* $*: Bootstrap simulations \
				 \n**********************************************\n"
	Rscript R/bootstrap.R $* $(num_cores)

output/%/blb.qs: R/blb.R output/%/setup.rda output/%/mle.qs
	@echo "\n**********************************************\
				 \n* $*: Bag of Little Bootstraps simulations \
				 \n**********************************************\n"
	Rscript R/blb.R $* $(num_cores)

output/%/armb.qs: R/armb.R output/%/setup.rda output/%/mle.qs
	@echo "\n**********************************************\
				 \n* $*: Averaged Robbins-Monro simulations \
				 \n**********************************************\n"
	Rscript R/armb.R $* $(num_cores)

sims_plots:
	@echo "\n*****************************************\
				 \n* Read simulation results and save plots \
				 \n******************************************\n"
	Rscript R/sims_plots.R

real_data: data/aps_dataset.zip  data/wilt_dataset.zip  output/wilt.pdf
#output/aps.pdf
data/aps_dataset.zip:
	@echo "\n**********************************************\
				 \n* APS SCANIA DATASET: Download \
				 \n**********************************************\n"
	mkdir -p data/aps/
	curl -o data/aps/aps_dataset.zip -L https://archive.ics.uci.edu/static/public/421/aps+failure+at+scania+trucks.zip
	unzip -o data/aps/aps_dataset.zip -d data/aps/

data/wilt_dataset.zip:
	@echo "\n**********************************************\
				 \n* WILT DATASET: Download \
				 \n**********************************************\n"
	mkdir -p data/
	curl -o data/wilt_dataset.zip -L https://archive.ics.uci.edu/static/public/285/wilt.zip
	unzip -o data/wilt_dataset.zip -d data/wilt/

output/aps.pdf: R/realdata/realdata_aps.R
	@echo "\n***********************************************\
				 \n* APS SCANIA DATASET: Run experiments \
				 \n***********************************************\n"
	mkdir -p output/
	Rscript R/realdata/realdata_aps.R

output/wilt.pdf: R/realdata/realdata_appendix.R
	@echo "\n**********************************************\
				 \n* WILT, MAGIC and LIKE DATASETS: Run experiment \
				 \n**********************************************\n"
	mkdir -p output/
	Rscript R/realdata/realdata_appendix.R