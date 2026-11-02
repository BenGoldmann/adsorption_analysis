###################
### DEFINITIONS ###
###################

# Interpreter for running the scripts
INTERPRETER = python3

# Directory containing the scripts and data files
HELPER_SCRIPTS_DIR = /home/bgoldman/Documents/mlcc_project/z_scripts/helper_scripts
ADSORPTION_SCRIPTS_DIR = /home/bgoldman/Documents/mlcc_project/z_scripts/adsorption_analysis

# Scripts
REMOVE_SOLVENT_FEEDER_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/remove_solvent_feeder.py
REMOVE_SOLVENT_SCRIPT = $(HELPER_SCRIPTS_DIR)/remove_solvent.py
ATTACHMENT_ANALYSIS_FEEDER_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/attachment_analysis_feeder.py
ATTACHMENT_ANALYSIS_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/attachment_analysis.py
STATISTICS_CALCULATOR_FEEDER_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/statistics_calculator_feeder.py
STATISTICS_CALCULATOR_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/statistics_calculator.py
DENSITY_CALCULATOR_FEEDER_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/density_calculator_feeder.py
DENSITY_CALCULATOR_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/density_calculator.py
SURFACE_SORTER_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/surface_sorter.py
ISOLATE_ATTACHED_FEEDER_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/isolate_attached_feeder.py
ISOLATE_ATTACHED_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/isolate_attached.py
PLOT_ATTACHMENT_VS_TIME_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/plot_attachment_vs_time.py
PLOT_DENSITY_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/plot_density.py
SNAPSHOT_CREATOR_FEEDER_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/snapshot_creator_feeder.py
SNAPSHOT_CREATOR_SCRIPT = $(ADSORPTION_SCRIPTS_DIR)/snapshot_creator.py

# Input and data files
INFO_FILE = SYSTEMINFO
ELEMENTS_DATA = $(ADSORPTION_SCRIPTS_DIR)/elements.data
DUMP_FILE = dump.lammpstrj

# Parameter files
REMOVE_SOLVENT_PARAMETERS_FILE = remove_solvent_parameters.in
ATTACHMENT_ANALYSIS_PARAMETERS_FILE = attachment_analysis_parameters.in
STATISTICS_CALCULATOR_PARAMETERS_FILE = statistics_calculator_parameters.in
DENSITY_CALCULATOR_PARAMETERS_FILE = density_calculator_parameters.in
ISOLATE_ATTACHED_PARAMETERS_FILE = isolate_attached_parameters.in
ISOLATE_ATTACHED_PARAMETERS_FILE = isolate_attached_parameters.in

# Output files
DUMP_FILE_NOSOLVENT = dump_nosolvent.lammpstrj
ATTACHMENT_MASTER = attachment_master.pkl
ATTACHMENT_SURFACE = attachment_surface.pkl
ATTACHMENT_STATISTICS_FILE = attachment_statistics.txt
ATTACHMENT_VS_TIME_FILE = attachment_vs_time.txt
TOTAL_DENSITY_FILE = total_density.txt
ATTACHED_DENSITY_FILE = attached_density.txt
ATTACHMENT_VS_TIME_PLOT = attachment_vs_time.jpg
SATURATED_SNAPSHOT_FILE = saturated_snapshot.jpg


#############
### RULES ###
#############

# Use the complete attachment-analysis chain as Make's default target.
# `:=` expands the variable immediately; `.DEFAULT_GOAL` selects `all` when
# make is run without an explicit target.
.DEFAULT_GOAL := all

# `wildcard` expands these shell-like patterns against directories under the
# current working directory. `sort` orders the matches and removes duplicates.
SIM_DIRS := $(sort $(wildcard \
	nio_*_dhta nio_*_terpineol \
	bao_*_dhta bao_*_terpineol \
	tio2_*_dhta tio2_*_terpineol \
	batio3_*_dhta batio3_*_terpineol))
# `foreach` checks each candidate and retains it only when its dump exists.
SIM_DIRS := $(foreach dir,$(SIM_DIRS),$(if $(wildcard $(dir)/$(DUMP_FILE)),$(dir)))

# `addsuffix` creates one target path per simulation directory.
REMOVE_SOLVENT_PARAMETER_FILES := $(addsuffix /$(REMOVE_SOLVENT_PARAMETERS_FILE),$(SIM_DIRS))
DUMP_FILES_NOSOLVENT := $(addsuffix /$(DUMP_FILE_NOSOLVENT),$(SIM_DIRS))
ATTACHMENT_ANALYSIS_PARAMETER_FILES := $(addsuffix /$(ATTACHMENT_ANALYSIS_PARAMETERS_FILE),$(SIM_DIRS))
ATTACHMENT_MASTER_FILES := $(addsuffix /$(ATTACHMENT_MASTER),$(SIM_DIRS))
ATTACHMENT_SURFACE_FILES := $(addsuffix /$(ATTACHMENT_SURFACE),$(SIM_DIRS))
STATISTICS_CALCULATOR_PARAMETER_FILES := $(addsuffix /$(STATISTICS_CALCULATOR_PARAMETERS_FILE),$(SIM_DIRS))
ATTACHMENT_STATISTICS_FILES := $(addsuffix /$(ATTACHMENT_STATISTICS_FILE),$(SIM_DIRS))
ATTACHMENT_VS_TIME_FILES := $(addsuffix /$(ATTACHMENT_VS_TIME_FILE),$(SIM_DIRS))
DENSITY_CALCULATOR_PARAMETER_FILES := $(addsuffix /$(DENSITY_CALCULATOR_PARAMETERS_FILE),$(SIM_DIRS))
TOTAL_DENSITY_FILES := $(addsuffix /$(TOTAL_DENSITY_FILE),$(SIM_DIRS))
ATTACHED_DENSITY_FILES := $(addsuffix /$(ATTACHED_DENSITY_FILE),$(SIM_DIRS))

# `firstword` and `subst` extract the leading surface prefix from each
# simulation-directory name. BaTiO3 produces both a BaO and a TiO2 output set;
# single-surface simulations produce only the matching prefix set.
surface_prefix = $(firstword $(subst _, ,$(1)))
sorted_outputs = $(if $(filter batio3_%,$(1)), \
	$(1)/tio2_attachment_vs_time.txt $(1)/bao_attachment_vs_time.txt \
	$(1)/tio2_total_density.txt $(1)/bao_total_density.txt \
	$(1)/tio2_attached_density.txt $(1)/bao_attached_density.txt, \
	$(1)/$(call surface_prefix,$(1))_attachment_vs_time.txt \
	$(1)/$(call surface_prefix,$(1))_total_density.txt \
	$(1)/$(call surface_prefix,$(1))_attached_density.txt)
SORTED_OUTPUT_FILES := $(foreach dir,$(SIM_DIRS),$(call sorted_outputs,$(dir)))
SORTED_ATTACHMENT_VS_TIME_FILES := $(foreach dir,$(SIM_DIRS),$(if $(filter batio3_%,$(dir)), $(dir)/tio2_attachment_vs_time.txt $(dir)/bao_attachment_vs_time.txt, $(dir)/$(call surface_prefix,$(dir))_attachment_vs_time.txt))
SORTED_TOTAL_DENSITY_FILES := $(foreach dir,$(SIM_DIRS),$(if $(filter batio3_%,$(dir)), $(dir)/tio2_total_density.txt $(dir)/bao_total_density.txt, $(dir)/$(call surface_prefix,$(dir))_total_density.txt))
SORTED_ATTACHED_DENSITY_FILES := $(foreach dir,$(SIM_DIRS),$(if $(filter batio3_%,$(dir)), $(dir)/tio2_attached_density.txt $(dir)/bao_attached_density.txt, $(dir)/$(call surface_prefix,$(dir))_attached_density.txt))

# plot_density.py scans below the SYSTEMINFO directory and creates one image for
# each solvent represented by these simulations. Its outputs depend on the
# available data, so declare only the solvent-specific images that can be made.
DENSITY_PLOT_FILES := $(strip \
	$(if $(filter %_dhta,$(SIM_DIRS)),$(CURDIR)/attached_density_dhta.jpg $(CURDIR)/total_density_dhta.jpg) \
	$(if $(filter %_terpineol,$(SIM_DIRS)),$(CURDIR)/attached_density_terpineol.jpg $(CURDIR)/total_density_terpineol.jpg))

# Saturated dump names depend on the simulation surface. BaTiO3 produces a
# separate top TiO2 snapshot and bottom BaO snapshot; symmetric slabs produce
# one snapshot with their own surface prefix.
saturated_dumps = $(if $(filter batio3_%,$(1)), \
	$(1)/tio2_saturated.lammpstrj $(1)/bao_saturated.lammpstrj, \
	$(1)/$(call surface_prefix,$(1))_saturated.lammpstrj)
SATURATED_DUMP_FILES := $(foreach dir,$(SIM_DIRS),$(call saturated_dumps,$(dir)))
SATURATED_SIDE_FILES := $(patsubst %.lammpstrj,%.side,$(SATURATED_DUMP_FILES))
ISOLATE_ATTACHED_PARAMETER_FILES := $(addsuffix /$(ISOLATE_ATTACHED_PARAMETERS_FILE),$(SIM_DIRS))
SNAPSHOT_PARAMETER_FILES := $(patsubst %_saturated.lammpstrj,%_snapshot_parameters.in,$(SATURATED_DUMP_FILES))
SNAPSHOT_IMAGE_FILES := $(patsubst %_saturated.lammpstrj,%_saturated_side.tga,$(SATURATED_DUMP_FILES)) $(patsubst %_saturated.lammpstrj,%_saturated_top.tga,$(SATURATED_DUMP_FILES))

# `define` stores a block of Make syntax as a reusable template. `$(1)` is the
# simulation-directory argument supplied later by `call`.
# Each rule has the form `target: prerequisites`; its tab-indented recipe runs
# only when the target is missing or older than one of its prerequisites.
# Concrete targets (rather than generic pattern-rule intermediates) are used so
# Make does not automatically remove these generated files after the build.
# In these recipes, `cd` lets each feeder infer the simulation directory.
# `$$$$` survives the two Make expansion passes in `eval`/`call` as `$`, allowing
# the shell to expand the existing `PYTHONPATH` environment variable. The prefix
# adds the project root (two levels above helper_scripts) to Python's import path,
# so `from z_scripts...` imports work even though the command runs inside a
# simulation subdirectory. `$(INTERPRETER)` only selects Python; it does not add
# project modules to Python's search path.
define SIMULATION_RULES
$(1)/$(REMOVE_SOLVENT_PARAMETERS_FILE): $(INFO_FILE) $(REMOVE_SOLVENT_FEEDER_SCRIPT) $(HELPER_SCRIPTS_DIR)/get_systeminfo_path.py $(HELPER_SCRIPTS_DIR)/read_parameters.py $(1)/$(DUMP_FILE)
	cd "$(1)" && PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$$${PYTHONPATH}" $(INTERPRETER) "$(REMOVE_SOLVENT_FEEDER_SCRIPT)"

$(1)/$(DUMP_FILE_NOSOLVENT): $(1)/$(REMOVE_SOLVENT_PARAMETERS_FILE) $(REMOVE_SOLVENT_SCRIPT)
	cd "$(1)" && $(INTERPRETER) "$(REMOVE_SOLVENT_SCRIPT)" "$(REMOVE_SOLVENT_PARAMETERS_FILE)"

$(1)/$(ATTACHMENT_ANALYSIS_PARAMETERS_FILE): $(INFO_FILE) $(ATTACHMENT_ANALYSIS_FEEDER_SCRIPT) $(HELPER_SCRIPTS_DIR)/get_systeminfo_path.py $(HELPER_SCRIPTS_DIR)/read_parameters.py $(1)/$(DUMP_FILE_NOSOLVENT)
	cd "$(1)" && PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$$${PYTHONPATH}" $(INTERPRETER) "$(ATTACHMENT_ANALYSIS_FEEDER_SCRIPT)"

# Grouped target syntax (`&:`) means this one recipe creates both pickle files.
# Make runs the analysis once if either output is missing or out of date, then
# treats both targets as updated; separate ordinary rules would run it twice.
$(1)/$(ATTACHMENT_MASTER) $(1)/$(ATTACHMENT_SURFACE) &: $(1)/$(ATTACHMENT_ANALYSIS_PARAMETERS_FILE) $(ATTACHMENT_ANALYSIS_SCRIPT)
	cd "$(1)" && PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$$${PYTHONPATH}" $(INTERPRETER) "$(ATTACHMENT_ANALYSIS_SCRIPT)" "$(ATTACHMENT_ANALYSIS_PARAMETERS_FILE)"

$(1)/$(STATISTICS_CALCULATOR_PARAMETERS_FILE): $(INFO_FILE) $(STATISTICS_CALCULATOR_FEEDER_SCRIPT) $(HELPER_SCRIPTS_DIR)/get_systeminfo_path.py $(HELPER_SCRIPTS_DIR)/read_parameters.py $(1)/$(ATTACHMENT_MASTER)
	cd "$(1)" && PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$$${PYTHONPATH}" $(INTERPRETER) "$(STATISTICS_CALCULATOR_FEEDER_SCRIPT)"

# This grouped target likewise records that statistics_calculator.py writes
# both the summary-statistics and per-frame time-series files in one run.
$(1)/$(ATTACHMENT_STATISTICS_FILE) $(1)/$(ATTACHMENT_VS_TIME_FILE) &: $(1)/$(STATISTICS_CALCULATOR_PARAMETERS_FILE) $(STATISTICS_CALCULATOR_SCRIPT) $(1)/$(ATTACHMENT_MASTER) $(1)/$(ATTACHMENT_SURFACE)
	cd "$(1)" && PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$$${PYTHONPATH}" $(INTERPRETER) "$(STATISTICS_CALCULATOR_SCRIPT)" "$(STATISTICS_CALCULATOR_PARAMETERS_FILE)"

$(1)/$(DENSITY_CALCULATOR_PARAMETERS_FILE): $(INFO_FILE) $(DENSITY_CALCULATOR_FEEDER_SCRIPT) $(HELPER_SCRIPTS_DIR)/get_systeminfo_path.py $(HELPER_SCRIPTS_DIR)/read_parameters.py $(1)/$(DUMP_FILE_NOSOLVENT)
	cd "$(1)" && PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$$${PYTHONPATH}" $(INTERPRETER) "$(DENSITY_CALCULATOR_FEEDER_SCRIPT)"

$(1)/$(TOTAL_DENSITY_FILE) $(1)/$(ATTACHED_DENSITY_FILE) &: $(1)/$(DENSITY_CALCULATOR_PARAMETERS_FILE) $(DENSITY_CALCULATOR_SCRIPT) $(1)/$(ATTACHMENT_MASTER) $(1)/$(ATTACHMENT_SURFACE) $(1)/$(ATTACHMENT_STATISTICS_FILE)
	cd "$(1)" && PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$$${PYTHONPATH}" $(INTERPRETER) "$(DENSITY_CALCULATOR_SCRIPT)" "$(DENSITY_CALCULATOR_PARAMETERS_FILE)"

# The sorter reads the three unsorted files from its working directory. Declare
# its exact prefix-specific outputs as a grouped target so one invocation
# creates them all (six for BaTiO3, three for the other surface types).
$(call sorted_outputs,$(1)) &: $(SURFACE_SORTER_SCRIPT) $(1)/$(ATTACHMENT_VS_TIME_FILE) $(1)/$(TOTAL_DENSITY_FILE) $(1)/$(ATTACHED_DENSITY_FILE)
	cd "$(1)" && PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$$${PYTHONPATH}" $(INTERPRETER) "$(SURFACE_SORTER_SCRIPT)"

$(1)/$(ISOLATE_ATTACHED_PARAMETERS_FILE): $(INFO_FILE) $(ISOLATE_ATTACHED_FEEDER_SCRIPT) $(HELPER_SCRIPTS_DIR)/get_systeminfo_path.py $(HELPER_SCRIPTS_DIR)/read_parameters.py $(1)/$(ATTACHMENT_VS_TIME_FILE) $(1)/$(ATTACHMENT_MASTER) $(1)/$(ATTACHMENT_SURFACE) $(1)/$(ATTACHMENT_STATISTICS_FILE)
	cd "$(1)" && PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$$${PYTHONPATH}" $(INTERPRETER) "$(ISOLATE_ATTACHED_FEEDER_SCRIPT)"

# The isolator selects the maximally attached frame(s) and writes matching
# `.side` metadata used later to orient top-selected snapshots for rendering.
$(call saturated_dumps,$(1)) $(patsubst %.lammpstrj,%.side,$(call saturated_dumps,$(1))) &: $(1)/$(ISOLATE_ATTACHED_PARAMETERS_FILE) $(ISOLATE_ATTACHED_SCRIPT) $(1)/$(DUMP_FILE_NOSOLVENT) $(1)/$(ATTACHMENT_VS_TIME_FILE) $(1)/$(ATTACHMENT_MASTER) $(1)/$(ATTACHMENT_SURFACE) $(1)/$(ATTACHMENT_STATISTICS_FILE)
	cd "$(1)" && PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$$${PYTHONPATH}" $(INTERPRETER) "$(ISOLATE_ATTACHED_SCRIPT)" "$(ISOLATE_ATTACHED_PARAMETERS_FILE)"

endef

# `call` fills in `$(1)` for each simulation; `eval` asks Make to parse the
# resulting text as real target rules.
$(foreach dir,$(SIM_DIRS),$(eval $(call SIMULATION_RULES,$(dir))))

# The snapshot feeder searches below this SYSTEMINFO directory, creates one
# parameter file beside each saturated dump, and the snapshot creator uses each
# file to render side/top VMD images. Each grouped rule runs its global script
# once even though it produces several per-simulation files.
ifneq ($(strip $(SNAPSHOT_PARAMETER_FILES)),)
$(SNAPSHOT_PARAMETER_FILES) &: $(SNAPSHOT_CREATOR_FEEDER_SCRIPT) $(INFO_FILE) $(HELPER_SCRIPTS_DIR)/get_systeminfo_path.py $(HELPER_SCRIPTS_DIR)/read_parameters.py $(SATURATED_DUMP_FILES) $(SATURATED_SIDE_FILES)
	PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$PYTHONPATH" $(INTERPRETER) "$(SNAPSHOT_CREATOR_FEEDER_SCRIPT)"

$(SNAPSHOT_IMAGE_FILES) &: $(SNAPSHOT_CREATOR_SCRIPT) $(SNAPSHOT_PARAMETER_FILES) $(SATURATED_DUMP_FILES) $(SATURATED_SIDE_FILES) $(ELEMENTS_DATA)
	PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$PYTHONPATH" $(INTERPRETER) "$(SNAPSHOT_CREATOR_SCRIPT)"
endif

# The plotter finds SYSTEMINFO from the Make working directory and recursively
# discovers the sorted files itself, so run it once globally (not once per
# simulation) and do not pass the raw, unsorted attachment time-series files.
$(ATTACHMENT_VS_TIME_PLOT): $(PLOT_ATTACHMENT_VS_TIME_SCRIPT) $(SORTED_ATTACHMENT_VS_TIME_FILES)
	PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$PYTHONPATH" $(INTERPRETER) "$(PLOT_ATTACHMENT_VS_TIME_SCRIPT)"

# This global plotter discovers sorted attached- and total-density files beneath
# the SYSTEMINFO directory. A grouped target reflects that one invocation writes
# all available solvent/profile-kind figures.
ifneq ($(DENSITY_PLOT_FILES),)
$(DENSITY_PLOT_FILES) &: $(PLOT_DENSITY_SCRIPT) $(SORTED_ATTACHED_DENSITY_FILES) $(SORTED_TOTAL_DENSITY_FILES)
	PYTHONPATH="$(abspath $(HELPER_SCRIPTS_DIR)/../..):$$PYTHONPATH" $(INTERPRETER) "$(PLOT_DENSITY_SCRIPT)"
endif

# Generate all.
# `.PHONY` marks `all` as a command-like goal rather than a file. Its prerequisites
# are all per-simulation attachment outputs, which trigger the rules above.
.PHONY : all
all : $(ATTACHMENT_MASTER_FILES) $(ATTACHMENT_SURFACE_FILES) $(ATTACHMENT_STATISTICS_FILES) $(ATTACHMENT_VS_TIME_FILES) $(TOTAL_DENSITY_FILES) $(ATTACHED_DENSITY_FILES) $(SORTED_OUTPUT_FILES) $(ATTACHMENT_VS_TIME_PLOT) $(DENSITY_PLOT_FILES) $(ISOLATE_ATTACHED_PARAMETER_FILES) $(SATURATED_DUMP_FILES) $(SATURATED_SIDE_FILES) $(SNAPSHOT_PARAMETER_FILES) $(SNAPSHOT_IMAGE_FILES)