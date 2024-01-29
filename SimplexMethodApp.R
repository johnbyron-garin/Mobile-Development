# John Byron Garin
# Exercise 10
# Simplex Method

library(shiny)

hasNegative <- function(initialTableau, nRow, nCol){
  hasNegative = FALSE
  n = 1
  while(n <= (nCol)){
    if(initialTableau[nRow, n] < 0){
      hasNegative = TRUE
      break
    }
    n = n + 1
  }
  return (hasNegative)
}

simplex <- function (initialTableau, isMax, problem){
  nRow = length(initialTableau[,1])
  nCol = length(initialTableau[1,])
  
  # identifying the columns in the right form to get their corresponding values
  
  # updating the basic solution based on the values per column
  if(isMax == TRUE){ # if the given problem is about maximization
    newColnames = initialTableau[1,]
    basicSolution = newColnames[-nCol]
    f = 1
    while(f < nCol){
      zeroCounter = 0
      oneCounter = 0
      for(elements in initialTableau[,f]){
        if(elements == 0){
          zeroCounter = zeroCounter + 1
        }else if(elements == 1){
          oneCounter = oneCounter + 1
        }else{
          basicSolution[f] = 0
          break
        }
      }
      if(zeroCounter == (nRow-1) & oneCounter == 1){
        rowIndex = 1
        while(rowIndex <= nRow){
          if(initialTableau[rowIndex,f] == 1){
            basicSolution[f]= initialTableau[rowIndex,nCol]
          }
          rowIndex = rowIndex + 1
        }
      }else{
        basicSolution[f] = 0
      }
      f = f + 1
    }
  }else{ # if the given problem is about minimization
    basicSolution = initialTableau[nRow,]
    basicSolution = basicSolution[-nCol]
    basicSolution[length(basicSolution)] = initialTableau[nRow,nCol]
  }
  
  optValue <- matrix(
    basicSolution,
    nrow = 1,
    ncol = (nCol-1),
    byrow = TRUE,
  )
  
  if(problem == TRUE){
    shippingNum <- matrix(
      optValue[1,9:23],
      nrow = 3,
      ncol = 5,
      byrow = TRUE,
      dimname = list(c("Denver", "Phoenix", "Dallas"), c("CA", "UT", "NM", "IL", "NYC"))
    )
    
    final = list(
      final.tableau = initialTableau,
      basic.solution = basicSolution,
      opt.val = optValue[1,(nCol-1)],
      shipping.num = shippingNum
    )
  }else{
    final = list(
      final.tableau = initialTableau,
      basic.solution = basicSolution,
      opt.val = optValue[1,(nCol-1)]
    )
  }
  
  return (final)
}

u <- shinyUI(pageWithSidebar(
  
  headerPanel("Simplex Method"),
  sidebarPanel(
    p(em("Vector")),
    textInput('matrix', 'Enter a vector to be placed by Row in a Matrix (comma delimited)', "3,2,1,0,0,0,66,9,4,0,1,0,0,180,2,10,0,0,1,0,200,-90,-75,0,0,0,1,0"),
    numericInput('nRow', 'Enter number of rows', 1, min = 1, max = 50),
    numericInput('nCol', 'Enter number of columns', 1, min = 1, max = 50),
    p(em("Row Names")),
    textInput('rowNames', 'Enter a vector of Row Names to be used in the Matrix (comma delimited)', "Production,Assembly,Quality Control,Unit Profit"),
    p(em("Column Names")),
    textInput('colNames', 'Enter a vector of Column Names to be used in the Matrix (comma delimited)', "x1,x2,s1,s2,s3,z,solution"),
    selectInput('var', 'Optimization:', c("Maximum" = TRUE, "Minimum" = FALSE)),
    actionButton("calc", "Calculate")
  ),
  
  mainPanel(
    h4('Answer'),
    verbatimTextOutput("simplex")
  )
))

s <- shinyServer(function(input, output) {
  
  calculate <- eventReactive(input$calc, {
    x <- c(as.numeric(unlist(strsplit(input$matrix,","))))
    rowNames <- c(unlist(strsplit(input$rowNames,",")))
    colNames <- c(unlist(strsplit(input$colNames,",")))
    nRow <- input$nRow
    nCol <- input$nCol
    isMax <- input$var
    problem <- FALSE
    if(!length(rowNames) == nRow | !length(colNames) == nCol){
      cat("ERROR!\nNumber of Row Names/Col Names is not equal to the number of Row/Column")
    }else{
      initialTableau <- matrix(
        x,
        nrow = nRow,
        ncol = nCol,
        byrow = TRUE,
        dimname = list(rowNames,colNames)
      )
      if(length(x) == nRow*nCol){
        cat("=====> Initial Tableau\n\n")
        print(initialTableau)
        
        #========================================================
        rows = length(initialTableau[,1])
        cols = length(initialTableau[1,])
        iteration = 1
        sw = 0
        while(hasNegative(initialTableau, rows, cols) != FALSE){
          # identifying the pivot column
          i = 1
          pivotColumnIndex = 1
          highestNegative = 0
          while(i <= (cols)){
            if(initialTableau[rows, i] < 0 & abs(initialTableau[rows, i]) > abs(highestNegative)){
              highestNegative = initialTableau[rows, i]
              pivotColumnIndex = i
            }
            i = i + 1
          }
          
          # computing the TR, and finding the lowest TR and the pivot row
          TR <- vector(mode="numeric", length=(rows-1))
          j = 1
          lowestTR = 10000000
          pivotRowIndex = 1
          while(j < rows){
            TR[j] = initialTableau[j,cols] / initialTableau[j,pivotColumnIndex]
            if(TR[j]==Inf || is.nan(TR[j])){
              TR[j] = 0
            }
            if(TR[j] > 0 & TR[j] < lowestTR){
              lowestTR = TR[j]
              pivotRowIndex = j
            }
            j = j + 1
          }
          
          # if the test ratio has no value greater than 0, then it will print no feasible solution
          if(all(TR<=0)) {
            sw = 1
            break
          }else{
            cat("\n=====> Iteration #")
            cat(iteration)
            cat("\n\n")
          }
          
          # computing for the normalized row and putting it into the tableau
          pivotElement = initialTableau[pivotRowIndex, pivotColumnIndex]
          normalizedRow = round((initialTableau[pivotRowIndex,] / pivotElement),4)
          initialTableau[pivotRowIndex,] = normalizedRow
          
          # updating the tableau with newRow after computing for it using the currentRow, normalizedRow, and c
          r = 1
          while(r <= (rows)){
            if(r != pivotRowIndex){
              currentRow = initialTableau[r,]
              c = initialTableau[r, pivotColumnIndex]
              newRow = currentRow - (normalizedRow * c)
              initialTableau[r,] = round(newRow,4)
            }
            r = r + 1
          }
          print(initialTableau)
          iteration = iteration + 1
        }
        if(sw == 0){
          simplex <- simplex(initialTableau, isMax, problem)
          cat("\n")
          
          cat("=====> Final Tableau\n\n")
          print(simplex[[1]])
          
          cat("\n")
          
          cat("=====> Basic Solution\n\n")
          print(simplex[[2]])
          
          cat("\n")
          
          cat("=====> Optimum\n\n")
          cat(simplex[[3]])
        }else{
          cat("\n\nIteration stops here,\nNo Feasible Solution!")
        }
      }else{
        cat("ERROR!\nNumber of input does not align with the number of rows and columns.\n\nSome elements may repeat/cut.\n\nNumber of elements in the matrix must be equal \nto the product of the number of rows and columns.")
      }
    }
    
  })
  
  output$simplex<-renderPrint({
    calculate()
  })
  
}
)
shinyApp(ui = u, server = s)