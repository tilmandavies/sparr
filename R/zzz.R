.onAttach <- function(...){
	s1 <- "   _____ ___  ____  ____  ____         \n"
	s2 <- "  / ___// _ \\/ _  \\/ __ \\/ __ \\        \n"
	s3 <- "  \\__ \\/ ___/ __  /  ___/  ___/        \n"
	s4 <- " ___/ / /  / / / / /\\ \\/ /\\ \\          \n"
	s5 <- paste0("/____/_/  /_/ /_/_/  \\__/  \\_\\   v", packageDescription("sparr")$Version, "\n\n")
	packageStartupMessage(paste("\n\nWelcome to\n",s1,s2,s3,s4,s5,"- type news(package=\"sparr\") for an overview\n- type help(\"sparr\") for documentation\n- type citation(\"sparr\") for how to cite\n",sep=""),appendLF=TRUE)
} #\n*type vignette(\"sparr2\") to access the accompanying article [unimplemented]\n*type citation(\"sparr2\") for how to cite use of this package\n"