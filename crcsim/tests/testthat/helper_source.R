################################################################################
# File: helper_source.R                                                        #
# Project: crc-sim                                                             #
# Created Date: 2026-05-05                                                     #
# Author: Matthew Carroll                                                      #
# -----                                                                        #
# Last Modified: 2026-05-15                                                    #
# Modified By: Matthew Carroll                                                 #
# -----                                                                        #
# Copyright (c) 2026 Syndemics Lab at Boston Medical Center                    #
################################################################################

# Prefer source-based loading during local test/debug runs so testthat::test_file()
# and debugger sessions exercise the current working tree, not an installed copy.
if (
	requireNamespace("pkgload", quietly = TRUE) &&
		file.exists("DESCRIPTION") &&
		dir.exists("R")
) {
	pkgload::load_all(
		path = ".",
		quiet = TRUE,
		export_all = FALSE,
		helpers = FALSE,
		attach_testthat = FALSE
	)
} else {
	library(crcsim)
}
