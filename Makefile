experiments = setting1 setting2

all: $(experiments)

clean:
	rm -r output

$(experiments): %: output/%/setup.rda output/%/gt.qs

output/%/setup.rda: R/%.R
	@echo "\n************************\n* Initialising $*: \n************************\n"
	mkdir -p output/$*/
	Rscript R/$*.R

output/%/gt.qs: R/ground_truth.R output/%/setup.rda
	@echo "\n*****************************\n* Computing the ground truth:\n*****************************\n"
	Rscript R/ground_truth.R $* 2

