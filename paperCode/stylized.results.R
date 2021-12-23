stylized.example.results <- lapply(list.files(path = "paperCode/supp_results"),
       function(x){readRDS(here::here("paperCode","supp_results",x))})

stylized.example.results <- do.call(rbind,stylized.example.results)
round(t(stylized.example.results),4)
