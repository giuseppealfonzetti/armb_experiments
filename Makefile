num_cores = 5

experiments = ber_id_fxd ber_eq_fxd ber_toe_fxd

all: $(experiments)

clean:
	rm -r output

$(experiments): %: output/%/setup.rda output/%/gt.qs output/%/mle.qs output/%/bt.qs output/%/blb.qs output/%/armb.qs

output/%/setup.rda: R/sett/%.R R/sett/common.R
	@echo "\n*********************************\
				 \n* $*: Initialisation\
				 \n*********************************\n"
	mkdir -p output/$*/
	Rscript R/sett/common.R
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






