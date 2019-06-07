## Special strains
WILDTYPES = config::get("constants")$wildtypes

## Feature categories
COMPUTED <- config::get("constants")$featureset$computed
DISCRETE <- config::get("constants")$featureset$discrete
REDUNDANT <- config::get("constants")$featureset$redundant

LEARNING_FEATURES <- config::get("constants")$featureset$learning

## Plotting constants
DEFAULT_PT_SIZE <- 6
DEFAULT_LN_WIDTH <- 0
DEFAULT_LN_SIZE <- 1
DEFAULT_AXIS_TEXT <- 16
DEFAULT_TITLE_TEXT <- 20
DEFAULT_TICK_TEXT <- 14
TITLE_SIZES <- theme(axis.title.x = element_text(size=DEFAULT_AXIS_TEXT),
                     axis.title.y = element_text(size=DEFAULT_AXIS_TEXT),
                     axis.text.x =  element_text(size=DEFAULT_AXIS_TEXT),
                     axis.text.y = element_text(size=DEFAULT_AXIS_TEXT),
                     title = element_text(size=DEFAULT_TITLE_TEXT),
                     legend.title = element_text(size=DEFAULT_AXIS_TEXT),
                     panel.border = element_rect(colour = "black", fill=NA, size=2))