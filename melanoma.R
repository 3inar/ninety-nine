library(magrittr)
library(dplyr)
library(stringr)

# Data downloaded from Meta-Kaggle
competitions <- readr::read_csv("data/meta_kaggle/competitions.csv")
competitions %<>% filter(str_detect(Slug, "melanoma",))

teams <- readr::read_csv("data/meta_kaggle/teams.csv")
teams %<>% filter(CompetitionId %in% competitions$Id) %>%
  filter(!is.na(PrivateLeaderboardSubmissionId))    # remove teams that entered but did not submit

submissions <- readr::read_csv("data/meta_kaggle/submissions.csv")
submissions %<>% filter(Id %in% teams$PrivateLeaderboardSubmissionId)

plot(PrivateScoreFullPrecision~PublicScoreFullPrecision,
     main="Melanoma Challenge:\nclear tendency for public over-optimism", 
     ylab="Private score", xlab="Public score",
     data=submissions, pch=".", xlim=c(.5,1), ylim=c(.5, 1))
abline(0,1)

