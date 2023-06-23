## ----Console.1, prompt = TRUE--------------------------------------------
10 + 5

(2 * ((10 + 5 - 3) * 2 - 9)) / 2^1

## ----Editor.1------------------------------------------------------------
15 + 5

((10 - 2)/4)^2

log(10) ## Takes the natural logarithm of the input.

log(exp(5)) ## Note the exponential function written as "exp()".

sin(pi/6) ## Sine function and π.

## ----Editor.2------------------------------------------------------------
# This is a command line which can be very long if you wish.
3 * 3 ## This is a code line with a command at the end.

# This is a comment line for a section.

## This is a comment for the subsection.

## ----Creating.Objects.1--------------------------------------------------
# Approach 1: try to use this approach.
x <- 6
print(x) ## Prints the object. I rarely use this function.
x ## Also prints the object.


z <- "R mini BootCamp" ## A character string. Note that strings in R are contained within double quotes.
z

r <- 1:10 ## ":" operator generates regular sequences.
r

# Approach 2: try to avoid this approach.
y = 4
y

w = "Go Wolfpack"
w

# Approach 3: try to avoid this approach.
assign("ncsu", "Go Wolfpack") ## Assings the "Go Wolfpack" character string to ncsu object.
ncsu

## ----Creating.Objects.2--------------------------------------------------
x <- 6 ## With lower case name.
x

X <- "NCSU" ## With capital letter name.
X

new.object <- 2:6 ## Don't use new_object or new-object.
new.object

conflicts(detail = TRUE)$.GlobalEnv

y <- "R mini BootCamp" ## Creating an object.
y
y <- 10:15 ## Overriding it with different values.
y

## ----R.Objects.Vectors.1-------------------------------------------------
# Creating an empty vector.
x <- vector("numeric", length = 10) ## Defines the class and length of a vector. You can use other classes that we will see in detail later.
x

# Vector with numeric value.
a <- c(5, 6, 7, 8, 9) ## Concatenate function which created a vector with numeric values.
a
length(a) ## Give the length of a vector.
dim(a) ## There is no dimension for vectors.
is.numeric(a) ## Checks whether the vector is numeric.
is.double(a) ## Numeric class in R is also called "double". So you can use "is.double" function as well.
class(a) ## Gives the class of an object in a character string.
str(a) ## Gives the details of object structure (class of the object and its values). Try to use it frequently, it is very useful.

seq(from = 1, to = 10, by = 1) ## seq function generates regular sequences.
seq(from = 1, to = 10, by = 4)
rep(1, each = 4) ## Replicates each value 3 times.
rep(c(2:5), times = 3) ## Replicates all value 4 times.
rep(c(2:5), each = 3) ## Replicates each value 3 times.

# Vector with integer values.
x <- c(1L, 2L, 3L, 4L) ## To create an integer in R use "L" after the numeric value.
x
is.integer(x) ## Checks whether the vector is integer.
class(x)
str(x)

# Vector with complex values.
a <- c(1 + 0i, 2 + 4i) ## Vector with complex values.
a
is.complex(a) ## Checks whether the vector is complex.
class(a)
str(a)

# Vector with character value.
a <- c("NCSU", "Wolfpack")
a
is.character(a) ## Checks whether the vector is character.
class(a)
str(a)

# Vector with logical values.
x <- c(TRUE, FALSE) ## Logical vector. See Logicals section for more details.
x
is.logical(x) ## Checks whether the vector is logical.
class(x)
str(x)

# Vector with names.
a <- c(1:3)
str(a)
attr(a, "names") <- c("First", "Second", "Third") ## A new attribute and description is added.
a ## Vector with names.
str(a) ## Named num.

b <- c(First = 1, Second = 2, Third = 3, Fourth = 4, Fifth = 5)
b ## Vector with names.
str(b) ## Named num.

## ----R.Objects.Vectors.2-------------------------------------------------
# Vector with numeric and character values.
a <- c(5, 6, 7, 8, "d") ## Note tha the last element in vector a is a character value.
str(a) ## Note tha class of vector a.

# Vector with logical and character values.
b <- c("a", TRUE) ## Character vector.
str(b)

# Vector with logical and numeric values.
x <- c(TRUE, 2) ## Numeric (TRUE will be converted into number 1). See Logicals section for more details.
str(x)
str(c(FALSE, 2)) ## Numeric (FALSE will be converted into number 0).

# Numeric vs. Integer values.
a <- c(1, 2)
is.integer(a)
is.numeric(a)

b <- c(1:2)
is.integer(b)
is.numeric(b)

## ----R.Objects.Vectors.3-------------------------------------------------
a <- c(5:15)
b <- c(10:20)
c <- c(25:30)
x <- c(b, a) ## Combining two vectors.
x 
y <- c(a, b, c) ## Combining three vectors.
y

## ----R.Objects.Vectors.4-------------------------------------------------
a <- c(1:10)
2 + a ## Vectorized operation.
1/a

## ----R.Objects.Vectors.5-------------------------------------------------
# Vectors with same length.
a <- seq(from = 1, to = 10, by = 2)
a
b <- seq(from = 1, to = 15, by = 3)
b
a/b

# Vectors with different length.
x <- seq(from = 1, to = 10, by = 2)
x
y <- seq(from = 1, to = 10, by = 3)
y

x/y ## Note the last value in this vector.

## ----R.Objects.Vectors.6-------------------------------------------------
a <- rnorm(n = 10000, mean = 0, sd = 1) ## Random number generator for the normal distribution with the specified mean and standard deviation. This is standard normal distribution.
head(x = a, n = 5) ## Prints the first 5 elements of a vector.
tail(x = a, n = 5) ## Prints the last 5 elements of a vector.
mean(a) ## Mean.
var(a) ## Variance.
sum((a - mean(a))^2) / (length(a) - 1) ## Same as above.
sd(a) ## Standard deviation.
sqrt(var(a)) ## Same as above. sqrt function is for square root.
min(a) ## Minimum value.
max(a) ## Maximum value.
sum(a) ## Total of all elements in a vector.

b <- c(seq(from = 10, to = 20, by = 4), seq(from = 10, to = 20, by = 2))
b
sqrt(b) ## Square root.
log(b)
sort(b, decreasing = FALSE, na.last = TRUE) ## Sorts the value of a vector alphabetically. na.last puts the missing values at the end. We will see the missing values later.
unique(b) ## Gives you the unique values in a vector.
sort(unique(b)) ## Unique values are sorted.

## ----R.Objects.Vectors.7-------------------------------------------------
a <- seq(from = 1, to = 200, by = 3)
a

## ----R.Objects.Factors.1-------------------------------------------------
a <- c("yes", "yes", "no", "yes", "no") ## A character vector.
str(a)

b <- factor(x = a) ## Creates the factor and gives you the "Levels".
b
is.factor(b) ## Checks whether the object is a factor.
class(b)
str(b) ## Note that the levels are automatically identified by alphabetic order.
attributes(b) ## Gives the object's attributes.
levels(b) ## Gives the levels.
table(b) ## Gives the number of levels by factors.
unclass(b) ## Shows the factors in numbers (no:1, yes:2).

b <- factor(x = a, levels = c("yes", "no"))
b ## Note that the levels are identified with levels argument.
table(b) ## Gives the number of levels by factors.
unclass(b) ## Shows the factors in numbers (no:1, yes:2).

attr(b, "levels") <- c("Aye", "Nay") ## Changin the level by using the attr function.
b

# Factor with names
a <- c(1:3)
attr(a, "names") <- c("First", "Second", "Third") ## A new attribute and description is added.
b <- factor(x = a)
b ## Note that the levels are identified with levels argument.

## ----R.Objects.Factors.2-------------------------------------------------
gl(n = 2, k = 4, labels = c("yes", "no")) ## Creates a factor object with 2 levels and 8 replications.

## ----R.Objects.Factors.3-------------------------------------------------
a <- c(1:15)
a.factor <- cut(x = a, breaks = c(min(a), 6, 12, max(a))) ## Note that first value is not included in the interval.
table(a.factor)
a.factor <- cut(x = a, breaks = c(min(a) - 1, 12, max(a))) ## Now all values are included.
table(a.factor)

a.factor <- cut(x = a, breaks = c(6, 12, max(a)), include.lowest = TRUE) ## If the minumum value is not specified, you need to use include.lowest argument.
table(a.factor)

a.factor <- cut(x = a, breaks = c(min(a) - 1, 6, 12, max(a)), labels = c("1st Group", "2nd Group", "3rd Group")) ## Note that we also defined the labels.
table(a.factor)
str(a.factor)

## ----R.Objects.Logicals.1------------------------------------------------
a <- TRUE
a
is.logical(a) ## Checks whether the vector is logical.
class(a)
str(a)

b <- FALSE
str(b)

x <- "TRUE"
str(x)

y <- c(TRUE, FALSE, FALSE)
str(y)

## ----R.Objects.Logicals.2------------------------------------------------
36 == 36 ## Checks equality of two numeric objects.
6 * 6 == 30 + 6 ## Checks the equality of two seperate calculations.
TRUE == FALSE ## Checks the equaility of two logical objects.
"NCSU" == "ncsu" ## Checks equality of two character objects.

36 != 36 ## Checks non-equaility.
36 != 6 ## Checks non-equaility.
TRUE != FALSE

1 < 0 ## Smaller than.
1 > 0 ## Bigger than.

1 <= 0 ## Smaller than and equal to.
1 >= 0 ## Bigger than and equal to.

-2:2 >= 0 ## Elementwise evaluation. Note that the shorter object is recycled fully over the longer one since the longer object length is a multiple of shorter object length.

(-2:2 >= 0) & (-2:2 <= 0) ## Elementwise evalutation. Note that the lengths of two vectors are same.

((-2:2) >= 0) && ((-2:2) <= 0) ## Only evaluates the first elements in each vector.

(-3:3 >= 0) & (-2:2 <= 0) ## Elementwise evaluation. Note that the shorter object cannot be recycled fully over the longer one since the longer object length is a multiple of shorter object length.
length(-3:3) %% length(-2:2) ## The reminder.

1:6 %in% 3:10 ## Elements of the left object is checked individually whether it matches with any of the right object elements.

## ----R.Objects.Logicals.3------------------------------------------------
a <- 0
class(a) ## Object "a" is a numeric object.
a == FALSE ## But it is still considered as FALSE since its value is "0".

b <- 1
class(b) ## Object "b" is a numeric object.
b == TRUE ## But it is still considered as TRUE since its value is "1".

x <- 2
class(x)
x == TRUE

## ----R.Objects.Logicals.4------------------------------------------------
a <- c(1:4, 2:5)
a == 3 ## Note tha TRUE values.

a <- c(5:15) ## First vector.
b <- c(10:20) ## Second vector.
x <- c(b, a) ## First and second vector are combined.
x

sort(unique(x)) ## Unique values are sorted.

sort(unique(x)) == sort(unique(c(a, b))) ## Checking if the sorting the unique values worked. Note the order of vectors (here it is c(a, b) but c vector is defined as c(b, a)).
sum(!(sort(unique(x)) == sort(unique(c(a, b))))) == 0 ## A quick control structure using the logical object and vectorized operations. So, there is no need to check it item by item.

## ----R.Objects.Logicals.5------------------------------------------------
x <- TRUE
if (x == TRUE) { ## We will see the details of IF statements later.
    print("My first IF statement.")
}

if (x) {
    print("IF statement result is TRUE, so it will be printed.")
}

if (!x) {
    print("IF statement result is FALSE, so it won't be printed.")
} else {
    print("IF statement result is TRUE, so it will be printed.")
}

## ----R.Objects.Matrices.1------------------------------------------------
matrix(data = 1:6, nrow = 2, ncol = 3, byrow = FALSE) ## Default is to fill the matrix by column.
matrix(data = 1:6, nrow = 2) ## You can define the length of one dimension.
matrix(data = 1:2, nrow = 2, ncol = 3) ## If the elements are not enought the data vector will be recycled to fill the whole matrix.
matrix(0, 2, 3) ## Creates a zero matrix. Note that value 0 recycles.

matrix(data = 1:6, nrow = 2, ncol = 3, byrow = TRUE) ## Filled by row.

a <- matrix(data = 1:6, nrow = 2, ncol = 3, dimnames = list(c("row1", "row2"), c("col1", "col2", "col3"))) ## We can also define the dimension names. Note that you need to use list() function. Wee will see the details of lists later.
a
is.matrix(a) ## Checks whether the object is a matrix.
class(a)
attributes(a) ## Gives the object's attributes.
str(a)

x <- c(1:6) ## You can also create matrix by first creating a vector then assigning its dimension.
dim(x) <- c(2, 3) ## Note that the matrix is filled by columns.
x

## ----R.Objects.Matrices.2------------------------------------------------
x <- matrix(data = 9:4, nrow = 3, ncol = 2, dimnames = list(c("row1", "row2", "row3"), c("col1", "col2")))
x

dim(x) ## Dimensions of a matrix.
dimnames(x) ## Dimension names.
nrow(x) ## Number of rows.
ncol(x) ## Number of columns.
rownames(x) ## Gives the row names.
colnames(x) ## Gives the column names.

rowSums(x) ## Row sums.
colSums(x) ## Column sums.
rowMeans(x) ## Row means.
colMeans(x) ## Column means.

diag(x = 2, nrow = 3, ncol = 3) ## Extract or replace the diagonal of a matrix, or construct a diagonal matrix. Here it creates an identity matrix.
diag(x = 3, nrow = 2, ncol = 2) ## Creates a diagonal matrix with values of 3.
diag(x = 3, nrow = 2, ncol = 3) ## One additional column.
diag(x = 3, nrow = 4, ncol = 3) ## One additional row.

y <- diag(x = 3, nrow = 5, ncol = 5)
diag(y) ## Extracts the diagonal.

z <- matrix(0 , nrow = 3, ncol = 3)
diag(z) <- 1:3 ## Assigns the given values to diagonal.
z

lower.tri(z, diag = FALSE) ## Gives the lower triangle of a matrix in logical values. You can use the result of this function to subset the lower triangle.
upper.tri(z, diag = FALSE) ## Gives the upper triangle of a matrix.
upper.tri(z, diag = TRUE) ## Diagonal is included.

## ----R.Objects.Matrices.3------------------------------------------------
x <- matrix(data = 1:6, nrow = 2, ncol = 3)
y <- matrix(data = 6:11, nrow = 3, ncol = 2)
z <- matrix(data = 10:15, nrow = 2, ncol = 3)
w <- matrix(data = 10:13, nrow = 2, ncol = 2)

t(x) ## Transpose.
det(w) ## Determinant.
solve(w) ## inverse.
eigen(w) ## Eigen values.

x + 1 ## Scalar summation.
x + z ## Matrix summation.
x / z ## Division by element.

2 * x ## Scalar multiplication.
c(2, 10) * x ## Scalar multiplication by row.
c(1, 2, 3) * y

x %% z ## Matrix multiplication

crossprod(x) ## Cross product.
kronecker(x, y) ## Kronecker product.
c(1:5) %o% c(1:5) ## Outer product
matrix(1:5, 5, 1) %*% matrix(1:5, 1, 5) ## Same result as above.

## ----R.Objects.Matrices.4------------------------------------------------
x <- matrix(data = 1:6, nrow = 2, ncol = 3)
y <- matrix(data = 6:11, nrow = 3, ncol = 2)
z <- matrix(data = 10:15, nrow = 2, ncol = 3)

rbind(x, z)
rbind(x, 1, 0) ## Recycling by columns.

cbind(x, z)
cbind(x, t(y))
cbind(x, 1, 0) ## Recyling by rows.

## ----R.Objects.Arrays.1--------------------------------------------------
array(data = 1:6, dim = c(2, 3), dimnames = NULL) ## Two-dimensional array.
matrix(data = 1:6, nrow = 2, ncol = 3, dimnames = NULL) ## Same as above.

x <- array(data = 1:12, dim = c(3, 2, 2)) ## 3-dimensional arrays.
x
dim(x)
is.array(x) ## Checks whether the object is an array.
class(x)
str(x)
attributes(x) ## Gives the object's attributes.

y <- array(data = 1:12, dim = c(3, 2, 2), dimnames = list(c("Row.1", "Row.2", "Row.3"), c("Col.1", "Col.2"), c("Dim.1", "Dim.2"))) ## Dimension names should be in a list form. See the Lists section first and check this function again.
y
str(y)

## ----R.Objects.Data.Frames.1---------------------------------------------
my.data <- data.frame() ## Creates an empty data frame.
is.data.frame(my.data) ## Checks whether the object is a data frame.
class(my.data)
str(my.data)

## ----R.Objects.Data.Frames.2---------------------------------------------
sample.size <- 30 ## Defining the sample size for later use.
column.1 <- round(rnorm(n = sample.size, mean = 5, sd = 1), digits = 2) ## Random number generator for the normal distribution with the specified mean and standard deviation. ## round function rounds the numeric values.
column.2 <- sample(x = c(-50:50, NA), size = sample.size, replace = TRUE, prob = NULL) ## Sample function takes a sample of the specified size from the elements of x using either with or without replacement. Creates an integer class.
column.3 <- sample(x = c("NCSU", "CALS", "Economics"), size = sample.size, replace = TRUE, prob = NULL)
column.4 <- factor(sample(x = c("Yes", "No"), size = sample.size, replace = TRUE, prob = NULL))
column.5 <- sample(x = c(TRUE, FALSE), size = sample.size, replace = TRUE, prob = NULL)

my.data <- data.frame(Column.1 = column.1, Column.2 = column.2, Column.3 = column.3, Column.4 = column.4, Column.5 = column.5) ## Creating data frame from scratch.
my.data
is.data.frame(my.data)
class(my.data)
str(my.data) ## Note that column.3 is factor variable but we wanted a character class.
attributes(my.data) ## Gives the object's attributes.

## ----R.Objects.Data.Frames.3---------------------------------------------
new.data <- data.frame(Column.1 = column.1, Column.2 = column.2, Column.3 = column.3, Column.4 = column.4, Column.5 = column.5, stringsAsFactors = FALSE) ## Creating data frame from scratch.
new.data
is.data.frame(new.data)
class(new.data)
str(new.data) ## Note that column.3 is character class which is what we wanted.

## ----R.Objects.Data.Frames.4---------------------------------------------
data.1 <- data.frame(c(1:5), c(6:10), c(11:15), c(16:20), stringsAsFactors = FALSE) ## Creating data frame from scratch without specific column names.
data.1

dim(data.1) ## Dimensions of a data frame.
dimnames(data.1) ## Dimension names.
nrow(data.1) ## Number of rows.
ncol(data.1) ## Number of columns.

rownames(data.1) ## Gives the row names.
colnames(data.1) ## Gives the column names.
names(data.1) ## Gives the column names.

column.names <- paste("Column", ".", 1:ncol(data.1), sep = "") ## Creating the generic column names automatically. paste function pastes the supplied values with a given string using vectorized operations.
column.names <- paste0("Column", ".", 1:ncol(data.1)) ## Same result as above. paste0 function pastes the supplied values with nothing using vectorized operations.
column.names
colnames(data.1) <- column.names ## Assignes the column names to the data frame by using the colnames.

# View(data.1) ## View the data frame in a new tab in interactive R session.
head(x = data.1, n = 2) ## Prints the first 2 elements of a data frame.
tail(x = data.1, n = 2) ## Prints the last 2 elements of a data frame.

rowSums(data.1) ## Row sums.
colSums(data.1) ## Column sums.
rowMeans(data.1) ## Row means.
colMeans(data.1) ## Column means.

## ----R.Objects.Lists.1---------------------------------------------------
# List without names.
a <- list(1, "a", TRUE, 1 + 4i)
a
is.list(a) ## Checks whether the object is a list.
class(a)
str(a)

# List with names.
a <- list(Numeric = 1, Character = "a", Logical = TRUE, Complex = 1 + 4i)
a
names(a) ## Gives a character vector of all the names of objects in a list.

# List inside of a list.
a <- list(c(2:4), "k", TRUE, list(rep(1, 3), rep(2, 4)))
a
is.list(a) ## Checks whether the object is a list.
str(a)

# List consists of different objects.
x <- c(1:2) ## Numeric vector.
y <- c("NCSU","Wolfpack", "Economics") ## Character vector.
z <- c(TRUE, FALSE, TRUE, FALSE, FALSE) ## Logical vector.
w <- factor(c("yes", "no", "no", "yes")) ## Factor vector.
v <- c(1 + 4i, 4 + 6i, 3 + 3i, 2 + 5i) ## Vector for complex values.
a <- matrix(data = 1:4, nrow = 2, ncol = 2, byrow = FALSE) ## Matrix.
b <- array(1:8, dim = c(2, 2, 2), dimnames = NULL) ## Array
data.1 <- data.frame(Column.1 = c(1:3), Column.2 = c(4:6), Column.3 = c(7:9), stringsAsFactors = FALSE)
my.list <- list(3, x, y, z, w, v, a, b, data.1) ## The list contains diffrent class of objects.
my.list
str(my.list)

## ----Working.Directory.1, eval = FALSE, tidy = FALSE---------------------
## # R code chunk is not evaluated.
## 
## R.home() ## Gives you the home directory of R software itself.
## 
## getwd() ## Gives the current working directory.
## my.current.dir <- getwd() ## Assigns the current working directory to an object.
## 
## setwd("Path of Working Directory") ## Sets the working directory to a new one. Note that this can be a relative path or a full path.
## setwd(my.current.dir) ## Using the assigned object, setting the working directory.
## 
## setwd("~") ## Changes the working directory to home directory.
## setwd("../") ## Double dots are used for moving up in the folder hierarchy.
## setwd("./") ## A single dot represents the current directory itself.
## setwd("/") ## Forward slash changes the working directory to the root.

## ----Workspace.1---------------------------------------------------------
x <- "NSCU"
class(x)  ## Gives the class of an object in a character string.
str(x) ## Gives the details of object structure (class of the object and its values). Try to use it, very useful.

attributes(x) ## This object does not have any attibutes yet.
attr(x, "Awesomeness Level") <- "Top Notch" ## A new attribute and description is added.
attributes(x) ## Not it has an user assign attribute.
structure(x, new.attribute = "This is a new attribute") ## Returns a new object with modified attributes.

attributes(cars) ## cars is a dataset from the datasets package in R.
class(cars)
str(cars)

## ----Workspace.2, eval = FALSE, tidy = FALSE-----------------------------
## # R code chunk is not evaluated.
## 
## ls() ## Shows the user created objects in your workspace.
## objects() ## Shows the user created objects in your workspace.
## ls.str() ## Shows the details of all objects in your workspace.
## 
## rm(x) ## Removes the object from your workspace.
## rm(c(x, y)) ## Removes multiple objects at the same time from your workspace.
## rm(list = ls()) ## Removes all objects from your workspace.

## ----Workspace.3, eval = FALSE-------------------------------------------
## # R code chunk is not evaluated.
## 
## help(options)
## options() ## View current option settings.
## options(digits = 3) ## Change an option setting. Number of digits to print on output.
## 
## history() # Displays last 25 commandss
## history(max.show = Inf) ## Displays all previous commandss
## 
## q() ## Ends R session. You will be prompted to save the workspace.

## ----File.System.1, eval = FALSE, tidy = FALSE---------------------------
## # R code chunk is not evaluated.
## 
## dir() ## Shows the files and folders in the working directory.
## dir.create("./folder") ## Will create a directory if it doesn't exist.
## file.exists("./RStudi_Setup.R") ## Will check to see if the directory exists.
## file.remove("./file.csv") ## Deletes the file or folder in the given path.
## unlink("./data.R") ## Deletes the file(s) or directories specified.
## 
## list.files("./") ## Lists the files in the given directory.
## list.files(pattern = "(.Rmd)") ## Lists the files with the selected pattern in the given directory.
## 
## if (!file.exists("./file.txt")) {
##     dir.create("./folder") ## Chekcs if the file exists, if not then creates the folder called "data".
## }
## file.remove("./file.xlsx") ## Removes the file or folder in the given path.
## 
## unzip("./file.zip") ## Extract the file from a zip archive.
## 
## ?files ## See the help file for more information such as renaming, appending and copying files.
## 
## fake.data <- rnorm(n = 1e6, mean = 0, sd = 1) ## Creates a fake data for size checking.
## object.size(fakeData) ## Gives the size of the R object in bytes.
## print(object.size(fakeData), units = "Mb") ## Gives the size of the R object in MB.
## file.info("../../R mini BootCamp.Rproj") ## Gives the file information of a file.
## file.info("../../R mini BootCamp.Rproj")$size ## Size of the file.

## ----File.System.2, eval = FALSE-----------------------------------------
## # R code chunk is not evaluated.
## 
## source("my_script.R") ## Loads "my_script.R" file from the current working directory.
## source("./my_script.R") ## Loads "my_script.R" file from the current working directory.
## source("../A Folder/my_script.R") ## Loads "my_script.R" file from another folder. See File System section for more information about file paths in R.
## source("FULL PATH of my_script.R file") ## Loads "my_script.R" with the full path.

## ----File.System.3, eval = FALSE-----------------------------------------
## # R code chunk is not evaluated.
## 
## x <- rbinom(n = 10, size = 1, prob = 0.5) ## Random number generator for the binomial distribution with parameters size and prob. n = number of observations, size = number of trials, prob = probability of success on each trial.
## y <- c("NCSU", "Wolfpack")
## save(x, file = "./x_object.RData") ## Saves the given objects to the RData file.
## save(x, y, file = "./xy_object.RData") ## You can save multiple objects in a Rdata file at the same time.
## rm(c(x, y)) ## Removes x and y objects.
## load("./xy_object.RData") ## Loads the RData file. Note that it overwrites the existing x and y objects if they are still in your workspace..
## 
## save.image() ## Saves the current workspace as .RData file. Note that it save as a hidden file.
## 
## saveRDS(object = x, file = "./xy_object.RData") ## Saves a single R object to a file.
## new.x <- readRDS(file = "xy_object.RData", refhook = NULL) ## Loads the x object as a new object.
## new.x <- load("./x.RData") ## Gives error.
## 
## dump(c("x", "y"), file = "./data.R") ## Dump can be used for multiple objects.
## rm(x)
## source("./data.R") ## Loads and runs the R code in the script.
## 
## dput(x, file = "./data.R") ## Dput can be used on single R objects.
## rm(x)
## new.x <- dget("./data.R") ## Loads and runs the R code in the script, then assigns a new name.
## 

## ----Getting.Help.1, eval = FALSE, tidy = FALSE--------------------------
## # R code chunk is not evaluated.
## 
## help(lm) ## Opens the help page for "lm" function which is for fitting linear models.
## ?lm
## ?"lm"
## ??lm ## Gives the search results for word "lm".
## ??errorsarlm ## If the package that contains the function is not installed, then you should use "??".
## ?":" ## Help for operator.
## ?"%in%"
## help.start() ## Opens the main page for R Help.
## help.search("covariance") ## Gives the search results for word "covariance".
## RSiteSearch("vecm") ## Opens your browser and searches for "vecm" on http://search.r-project.org.
## 
## find("lm") ## Tells you what package the function is in.
## apropos("lm") ## Returns a character vector giving the names of all objects in the search list that match your query.
## 
## args(lm) ## Presents the arguments of the function.
## example(lm) ## Presents an example of the searched function.
## demo(graphics) ## Gives a user-friendly interface to run some demonstration R scripts.

## ----Packages.1, eval = FALSE--------------------------------------------
## # R code chunk is not evaluated.
## 
## install.packages("tidyr") ## Installs single package.
## library("tidyr") ## Loads single package.
## 
## install.packages(c("RColorBrewer", "stringr")) ## Installs multiple packages.
## lapply(c("RColorBrewer", "stringr"), library, character.only = TRUE) ## Loads multiple packages.
## 
## packageVersion("tidyr") ## Current Version of the package.
## detach("package:RColorBrewer", unload = TRUE) ## Unloads package.

## ----Packages.2, eval = FALSE--------------------------------------------
## ## R chunk is not evaluated.
## 
## # Check the CRAN mirror.
## getOption("repos")
## 
## # Lists loaded packages in your global environment.
## (.packages())
## 
## # Some packages need to be installed from the source.
## install.packages("rgdal", type = "source")
## 
## # Installing a package from GitHub repository.
## install.packages("devtools") ## devtools package is necessary to install packages from GitHub repositories.
## devtools::install_github("tidyverse/ggplot2") ## user.name/package.name
## library("ggplot2")
## # devtools::install_github("hadley/devtools")
## 
## # Installing a package from bioconductor website.
## source("http://bioconductor.org/biocLite.R")
## biocLite("rhdf5")
## library("rhdf5")
## 
## # Installing CRAN Task Views.
## ## CRAN Task View (https://cran.r-project.org/web/views/) gives you the collection of packges in terms of area. For example spatial, econometrics, graphics and etc. To automatically install these views, the "ctv"" package needs to be installed, and then the views can be installed via "install.views" or "update.views" (which first assesses which of the packages are already installed and up-to-date) functions.
## install.packages("ctv") ## Intalling the necessary package to install views.
## library("ctv") ## Loading the package.
## 
## available.views(repos = NULL) ## Gives you all the views available.
## install.views("Econometrics") ## Installing the "Econometrics" view.
## update.views("Econometrics") ## Updating the "Econometrics" view.
## 
## # Updating packages from the editor or console.
## update.packages() ## Updates all packages from CRAN.
## 
## devtools::install_github("hrbrmstr/dtupdate") ## Updates Git sourced package. dtupdate package is used for this purpose.
## library("dtupdate")
## github_update() ## See what packages are avilable to update.
## 
## # Some other tools with packages.
## find.package("devtools") ## Shows you where (file loaction) the packge is installed in your computer.
## 
## search() ## Displays all the packages in the global environment.
## utils::installed.packages() ## Displays all the packages that are installed in your computer.
## available.packages() ## Displays all the R packages that are available.
## head(rownames(available.packages()), 3) ## Shows the names of the first 3 packages which are available.
## 
## # Package help.
## install.packages("sp") ## "sp" package is for spatial analysis.
## library("sp")
## vignette("sp") ## Opens the vignette for selected package if available.

## ----Packages.3, eval = FALSE--------------------------------------------
## # R code chunk is not evaluated.
## 
## Load.Install <- function(package_names) {
##     is_installed <- function(mypkg) is.element(mypkg, utils::installed.packages()[ ,1])
##     for (package_name in package_names) {
##         if (!is_installed(package_name)) {
##             utils::install.packages(package_name, dependencies = TRUE)
##         }
##         suppressMessages(library(package_name, character.only = TRUE, quietly = TRUE, verbose = FALSE))
##     }
## }
## 
## Load.Install(c("plyr", "dplyr", "tidyr", "sp"))

## ----Citation.1----------------------------------------------------------
# Cite R software.
citation()

# Cite R packages.
citation("ggplot2")

## ----Arithmatic.Operators, echo = FALSE----------------------------------
data.name <- "Arithmetic Operations"
data.colnames <- c("Operator", "Description")

data.list <- list(c("+", "Addition"),
                  c("-", "Substraction"),
                  c("*", "Multiplication"),
                  c("/", "Division"),
                  c("^ or **", "Exponentiation"),
                  c("%%", "Reminder"),
                  c("%/%", "Quotient"))

data <- as.data.frame(matrix(data = as.character(NA), nrow = length(data.list), ncol = length(data.colnames), byrow = FALSE), stringsAsFactors = FALSE)
colnames(data) <- data.colnames
for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
        data[i, j] <- data.list[[i]][j]
    }
}
knitr::kable(data, align = "c", caption = data.name)

## ----Logical.Operators, echo = FALSE-------------------------------------
data.name <- "Logical Operators"
data.colnames <- c("Operator", "Description")

data.list <- list(c("<", "Less than"),
                  c("<=", "Less than or equal to"),
                  c(">", "Greater than"),
                  c(">=", "Greater than or equal to"),
                  c("==", "Exactly equal to"),
                  c("!=", "Not equal to"),
                  c("| and ||", "OR"),
                  c("& and &&", "AND"),
                  c("%in%", "QLeft to rigth matching"))

data <- as.data.frame(matrix(data = as.character(NA), nrow = length(data.list), ncol = length(data.colnames), byrow = FALSE), stringsAsFactors = FALSE)
colnames(data) <- data.colnames
for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
        data[i, j] <- data.list[[i]][j]
    }
}
knitr::kable(data, align = "c", caption = data.name)

## ----Other.Operators, echo = FALSE---------------------------------------
data.name <- "Other Operators"
data.colnames <- c("Operator", "Description")

data.list <- list(c(":", "Generates regular sequences"))

data <- as.data.frame(matrix(data = as.character(NA), nrow = length(data.list), ncol = length(data.colnames), byrow = FALSE), stringsAsFactors = FALSE)
colnames(data) <- data.colnames
for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
        data[i, j] <- data.list[[i]][j]
    }
}
knitr::kable(data, align = "c", caption = data.name)

## ----Missing.Values.Vectors.1--------------------------------------------
x <- c(1:3, NA, 5:7, NA) ## Numeric vector.
x
str(x) ## Fourth value is missing.
is.na(x) ## Checks the missing values.
complete.cases(x) ## Checks the non-missing values.
is.na(x) == !complete.cases(x) ## Same functions.
!is.na(x) == complete.cases(x) ## Same functions.
sum(is.na(x)) ## Gives the number of missing values.
any(is.na(x)) ## Is there a missing value?

y <- c("a", "b", "c", NA, "NA") ## Character vector.
y
str(y)
is.na(y) ## Note that "NA" is not missing, it is a character string with values "NA".

## ----Missing.Values.Factors.1--------------------------------------------
a <- factor(x = c("yes", NA, "no", "yes", NA)) ## Creates the factor and gives you the 'Levels'.
a
str(a)
is.na(a)

## ----Missing.Values.Logicals.1-------------------------------------------
a <- c(TRUE, NA, FALSE, NA)
a
str(a)
is.na(a)

## ----Missing.Values.Matrices.1-------------------------------------------
a <- rep(c(3, NA, 2), each = 2)
b <- matrix(data = a, nrow = 2, ncol = 3)
b
str(b)
is.na(b)

na.check <- is.na(b)
class(na.check)
str(na.check)

## ----Missing.Values.Arrays.1---------------------------------------------
a <- sample(x = c(1:3, NA), size = 12, replace = TRUE, prob = NULL) ## Sample function takes a sample of the specified size from the elements of x using either with or without replacement. Creates an integer class.
b <- array(data = a, dim = c(2, 3, 2)) ## Two-dimensional array.
b
str(b)
is.na(b)

na.check <- is.na(b)
class(na.check)
str(na.check)

## ----Missing.Values.Data.Frames.1----------------------------------------
x <- data.frame(c(NA, 1:3, NA), c(NA, 4, NA, 5:6), c(7:9, NA, NA), c(10:14), stringsAsFactors = FALSE)  ## Creating data frame from scratch without specific column names.
colnames(x) <- paste0("Column", ".", 1:ncol(x)) ## Assignes the column names to the data frame by using the colnames.
x
str(x)
is.na(x)

sum(is.na(x)) ## Gives number of the missing values
any(is.na(x)) ## Is there a missing value?
colSums(is.na(x)) ## Missing values by columns.
rowSums(is.na(x)) ## Missing values by rows.

## ----Missing.Values.Lists.1----------------------------------------------
a <- list(NA, c(2:4), "k", rep(NA, 2), list(rep(NA, 3)))
a
is.na(a) ## Returns only 3 logical results. This is because is.na function thinks c(2:4) and rep(NA, 2) are single elements.
str(a)

b <- unlist(a)
b
str(b)
is.na(b) ## Unlisting gives the correct result.

## ----Missing.Values.Miscellaneous.1--------------------------------------
var(10) ## Variance of a number which returns NA.
sd(8) ## Standard deviation of a number which returns NA.
c(1, NA) == NA

NA + 1
NA * 2

x <- c(4, 5, NA)
x < 10
x + 2
x * 2
sum(x)
sum(x, na.rm = TRUE) ## Only the missing values are removed.
mean(x)
mean(x, na.rm = TRUE)
var(x)
var(x, na.rm = TRUE)

a <- rep(c(3, NA, 2), each = 2)
b <- matrix(data = a, nrow = 2, ncol = 3)
b
sum(b, na.rm = TRUE) ## Only the missing values are removed.
sum(b, na.omit = TRUE) ## The row of the missing value is removed.

length(x)

## ----Missing.Values.Miscellaneous.2--------------------------------------
x <- c(NA, NaN) ## Numeric vector.
str(x)
is.na(x) ## Both NA and NaN are considered as NA.
is.nan(x) ## Only NA is not considered as NaN.
is.nan(0/0)

## ----Subsetting.Vectors.1------------------------------------------------
# Vector with numeric values.
a <- c(1:10)
a

a[3] ## Selecting the 3rd element. Preserving subsetting with unnamed vectors.
a[[3]] ## Use it for vectors with names. Simplifying subsetting with unnamed vectors. Same result as above.

a[c(2, 3, 5)] ## Selecting the 2nd, 3rd and the 5th elements.
a[1:3]
a[10:5]
a[3:length(a)]
a[c(seq(1, 10, 2))]
a[c(2, 3, 5, 6)][c(1, 2)] ## Subsetting two times.
a[c(2, 3, 5, 6)][c(1, 2)][1] ## Subsetting three times..

# Negative indexing
a <- c(1:10)
a[-1]
a[-c(1:4)]
a[-c(1, 4)]

# Vector with names.
b <- c(First = 1, Second = 2, Third = 3, Fourth = 4, Fifth = 5)
b
b[1] ## Preseving subsetting with named vectors.
b[[1]] ## Simplifying subsetting with named vectors. Compared to above result, they are not the same.

b["First"]
b[["First"]] ## Simplifying subsetting.
b[[1]] ## Same as above.
b[c("First", "Third")]
b[c(1, 3)] ## Same as above.
# b[[c("First", "Third")]] ## Error. You cannot use multiple indices with "[[ ]]".Instead use the below command.
unname(b[c("First", "Third")]) ## You can use "unname" functio to drop the names.
c(b[["First"]], b[["Third"]])
b[c(names(b)[c(3, 5)])]
b[c(names(b)[-c(3, 5)])]

## ----Subsetting.Factors.1------------------------------------------------
a <- factor(x = c("yes", NA, "no", "yes", NA))  ## Creates the factor and gives you the 'Levels'.
a
a[3] ## Selecting the 3rd element. Preserving subsetting with unnamed factors.
a[[3]] ## Use it for factors with names. Simplifying subsetting with unnamed vectors. Same result as above.

a[c(1:3)]
a[1, drop = TRUE] ## Drops the unused levels.

a <- factor(x = c(First = "yes", Second = NA, Third = "no", Fourth = "yes", Fifth = NA))
a
a[1] ## Preseving subsetting with named factors.
a[[1]] ## Simplifying subsetting with named factors. Compared to above result, they are not the same.
unname(a[1]) ## Same as above.
a["First"] ## You can also use names.

a[1, drop = TRUE] ## Preserving subsetting. Drops the unused levels.
a[[1]][ , drop = TRUE] ## Simplifying subsetting. Drops the unused levels.

## ----Subsetting.Logicals.1-----------------------------------------------
y <- sample(x = c(TRUE, FALSE, NA), size = 5, replace = TRUE)
y
y[1]
y[[1]]
y[c(2:4)]

# Named logicals are dropped since they are rarely used.

## ----subsetting.matrices.1-----------------------------------------------
a <- matrix(data = 1:9, nrow = 3, ncol = 3, dimnames = list(c("row1", "row2", "row3"), c("col1", "col2", "col3")))
a

a[2] ## Gives the second element in a matrix. Note that to index of elements are by columns. Second element is on the first column second row.
a[[2]] ## Same result.

a[1, 1] ## Simplifying subsetting. First row and first column.
a[1, 1, drop = FALSE] ## Preserving subsetting.
a[1:2, 1:2]
a["row1", "col1"] ## You can also use names with quotes.

a[ , 1] ## Simplifying subsetting. First column.
a[ , 1, drop = FALSE] ## Preserving subsetting.
unname(a[, 1])
a[, 2:3]

a[1, ] ## Simplifying subsetting. First row.
a[1, , drop = FALSE] ## Preserving subsetting.
unname(a[1, ])
a[2:3, ]

a[-1, ] ## Negative subsetting
a[, -1]
a[-c(2:3), -1]

## ----Subsetting.Arrays.1-------------------------------------------------
a <- array(data = 1:12, dim = c(3, 2, 2)) ## 3-dimensional array.
a

a[12] ## Gives the 12th element.
a[[12]] ## Same result.

# a[1, 1] ## Incorrect number of dimensions
a[1, 1, 1] ## Simplifying subsetting.
a[1, 1, 1, drop = FALSE] ## Preserving subsetting.
a[1:2, 1:2, 1:2]

a[1, , ] ## Simplifying subsetting.
a[1, , , drop = FALSE] ## Preserving subsetting.
a[1:2, , ]

a[, 1, ] ## Simplifying subsetting.
a[, 1, , drop = FALSE] ## Preserving subsetting.
a[, 1:2, ]

a[ , , 1] ## Simplifying subsetting.
a[ , , 1, drop = FALSE] ## Preserving subsetting.
a[, , 1]

a[-1, , ] ## Negative subsetting
a[, -1, ]
a[-c(2:3), , -1]

## ----Subsetting.Data.Frames.1--------------------------------------------
x <- data.frame(c(1:5), c(6:10), c(11:15), c(16:20), stringsAsFactors = FALSE)  ## Creating data frame from scratch without specific column names.
colnames(x) <- paste0("Column", ".", 1:ncol(x)) ## Assignes the column names to the data frame by using the colnames.
x

# Subsetting columns.
x[, 1] ## Simplifying subsetting for columns. Column 1 values only.
x[[1]] ## Same as above. Note that it subsets the columns only.
x[["Column.1"]]
x[, "Column.1"]
x["Column.1"]
x$Column.1 ## Same as above. 
x$"Column.1" ## Same as above.
x[2:4, 1]

x[, 1, drop = FALSE] ## Preserving subsetting for columns. Column 1 values only.
x[1] ## Note that it subsets the columns only.
x[2:4, 1, drop = FALSE]

# Subsetting rows.
x[1, ] ## Structure of the data is preserved
unname(as.matrix(x[1, ])[1, ]) ## This is the simplified subsetting for rows. Note that as.matrix function coerce the subsetted data frame into matrix. We will see the details of coercion later.

# Subsetting row and columns.
x[1, 1]
x[2:4, 1]
x[c(1, 3), c(2, 4)]

# Negative subsetting
x[-1, -1] ## First row and first column is deleted.
x[-c(2:4), ] ## Row 2, 3, 4 are deleted.

## ----Subsetting.Lists.1--------------------------------------------------
a <- list(c(1:5), c("a", "b"), c(TRUE, FALSE), list(c(6:10), c("c, d")))
a

a[[1]] ## Simplifying subsetting.
a[1] ## Preserving subsetting. Note that the result is still a list.

a[[2]][1] ## [[]] helps us to get in the second element in the list. [] helps us the subset the second element in the list.

a[[4]][[1]]
a[[4]][[1]][2]

a <- list(Numeric = c(1:3), Character = c("a", "b"), Logical = c(TRUE, FALSE))
a
a[["Numeric"]]
a$Numeric
a$"Numeric"
a$Numeric[2]

## ----Conditional.Subsetting.1--------------------------------------------
# Vectors
x <- c(1:10)
x
a <- x > 4 ## This is our condition.
x[a]
x[!a]
x[x < 2 | x > 8]

# Matrices
a <- matrix(data = 1:9, nrow = 3, ncol = 3, dimnames = list(c("row1", "row2", "row3"), c("col1", "col2", "col3")))
a
a > 4 ## Condition.
a[a > 4] ## Condition is applied to all matrix elements.
b <- unname(a[, 2, drop = TRUE]) ## Second column.
b
b[b > 4] ## Condition on the second column.

# Data frames
x <- data.frame(c(1:5), c(6:10), c(11:15), c(16:20), stringsAsFactors = FALSE)  ## Creating data frame from scratch without specific column names.
colnames(x) <- paste0("Column", ".", 1:ncol(x)) ## Assignes the column names to the data frame by using the colnames.
x
x[x > 4] ## Condition is applied to all data frame elements.
x[x$Column.1 > 2, ]
x[x$Column.1 > 2 & x$Column.3 > 13, ]
x[x$Column.1 > 2 & x$Column.3 > 13, ]$Column.4
x[x$Column.1 > 2 & x$Column.3 > 13, c("Column.2", "Column.3")]

## ----Missing.Values.Problem.1--------------------------------------------
# Vectors
x <- sample(x = c(1:10, rep(NA, 5)), size = 10, replace = TRUE, prob = NULL)
x
a <- x > 4 ## When there is NA the condition produces NA.
a
b <- x[!is.na(x)] ## NA values are excluded.
b
b[b > 4] ## Values larger than 4.
x[which(x > 4)] ## You can use this directly. which functions omits the NA value automatically.

# Matrices
x <- sample(x = c(1:10, rep(NA, 5)), size = 9, replace = TRUE, prob = NULL)
a <- matrix(data = x, nrow = 3, ncol = 3, dimnames = list(c("row1", "row2", "row3"), c("col1", "col2", "col3")))
a
a > 4 ## When there is NA the condition produces NA.
a[which(a > 4)] ## Condition is applied to all matrix elements.
b <- unname(a[, 2, drop = TRUE]) ## Second column.
b
b[which(b > 4)] ## Condition on the second column.

# Data frames
x <- data.frame(c(NA, 1:3, NA), c(NA, 4, 10, 5:6), c(7:9, NA, NA), c(10:14), stringsAsFactors = FALSE)  ## Creating data frame from scratch without specific column names.
colnames(x) <- paste0("Column", ".", 1:ncol(x)) ## Assignes the column names to the data frame by using the colnames.
x
x > 4 ## When there is NA the condition produces NA.
is.na(x)
complete.cases(x) ## Gives the row with all non-NA values.
a <- x[complete.cases(x), ] ## Rows with non missing elements.
a ## Non-NA data frame.
a[a$Column.1 > 1, ] ## Apply the condition.
a[a$Column.1 > 1, c("Column.2", "Column.4")]

x[which(x$Column.1 > 1), ]
x[which(x$Column.1 > 1 & x$Column.3 > 7), ]
x[which(x$Column.1 > 2 & x$Column.4 > 11), ]$Column.4
x[which(x$Column.1 > 1 & x$Column.3 > 8), c("Column.2", "Column.3")]

## ----Assignment.by.Subsetting.1------------------------------------------
# Vectors
x <- c(1:10)
x
x[x > 4] <- NA
x
x[is.na(x)] <- 0

# Matrices
a <- matrix(data = 1:9, nrow = 3, ncol = 3, dimnames = list(c("row1", "row2", "row3"), c("col1", "col2", "col3")))
a
a[1, c(1:3)] <- NA
a

# Data frames
x <- data.frame(c(1:5), c(6:10), c(11:15), c(16:20), stringsAsFactors = FALSE)  ## Creating data frame from scratch without specific column names.
colnames(x) <- paste0("Column", ".", 1:ncol(x)) ## Assignes the column names to the data frame by using the colnames.
x
x[3, 4] <- 10000
x
x[x$Column.1 > 2, ] <- NA
x

## ----Coercion.1----------------------------------------------------------
x <- c(0:6)
class(x) ## The class of x is integer.

as.numeric(x) ## Coerces x as a numeric.
as.character(x) ## Coerces x as a character.
as.complex(x) ## Coerces x as a complex.
as.factor(x) ## Coerces x as a factor.
as.logical(x) ## Coerces x as a logical (0 is FALSE and everything greater than 0 is TRUE).
as.matrix(x) ## Coerces x as a matrix.
as.array(x) ## Coerces x as an arrray.
as.data.frame(x) ## Coerces x as a data frame.
as.list(x) ## Coerces x as a list.

## ----Dates.Times.1-------------------------------------------------------
d1 <- base::date() ## Current date and time.
d1
class(d1) ## Class is character

d2 <- Sys.Date() ## System date.
d2
class(d2) ## Class is "Date".

d3 <- Sys.time() ## System time.
d3
str(x) ## Class is "POSIXct"

## ----Dates.Times.2-------------------------------------------------------
x <- Sys.time()
x
str(x) ## Class is "POSIXct"

a <- as.POSIXlt(x) ## Coerced to as.POSIXlt.
str(a) ## Class is "POSIXltt"
names(unclass(a)) ## Gives the names, after it is unclassed.
a$sec
a$yday
a$mday

as.POSIXct(a) ## Coerced to as.POSIXct.

## ----Dates.Times.3-------------------------------------------------------
x <- as.Date("1970-01-01")
str(x) ## Class is "Date".

unclass(x) ## Note that the starting date is 1970-01-01 for date class.
unclass(as.Date("1970-01-02"))
unclass(as.Date("1969-12-31")) ## Dates before the starting date are represented by negative numbers.
unclass(as.Date(Sys.Date())) ## Number of days since the first date.

as.Date(0) ## The first date.
as.Date(c(-2:2))

## ----Dates.Times.4-------------------------------------------------------
x <- Sys.Date() ## System date.

weekdays(x, abbreviate = FALSE) ## The weekday.
months(d2, abbreviate = FALSE) ## The month.
julian(x) ## Gives the number of the days since the origin and the origin is given in the result.

format(x, "%a %b %d") ## Formats the date class object with the desired representation.
strftime(x, origin = "1970-01-01", tz = "UTC", format = "%B %d, %Y %H:%M") ## Similar to above.
strftime(x, origin = "1970-01-01", tz = "UTC", format = "%A %B %Y") ## Similar to above.
format(as.POSIXct("Feb 03, 2017 09:12 PM", tz = "UTC", format = "%b %d, %Y %I:%M %p"), "%Y%m%d%H%M")

date.string <- c("January 10, 2012 10:40", "December 9, 2011 09:10:00") ## A date string with a specific representation.
a <- strptime(date.string, format = "%B %d, %Y %H:%M") ## Note that the format needs match your date string format.
a
str(a)
as.Date(a)
# ?strptime ## Check the arg of strptime function.

## ----Dates.Times.5-------------------------------------------------------
year <- "2017"
d <- as.Date(paste0(year, "-01-01"))
tuesdays <- d + seq(by = 7, (2 - as.POSIXlt(d)$wday) %% 7, 364 + (months(d + 30 + 29) == "February")) ## Note that the day number starts on sunday with 0 and ends on saturday with 6.
tuesdays[tapply(seq_along(tuesdays), as.POSIXlt(tuesdays)$mon, max)]

d <- as.Date(paste0(year, "-01-01"))
saturdays <- d + seq(by = 7, (6 - as.POSIXlt(d)$wday) %% 7, 364 + (months(d + 30 + 29) == "February"))
saturdays[tapply(seq_along(saturdays), as.POSIXlt(saturdays)$mon, max)]

## ----Function.Call.1-----------------------------------------------------
log(10) ## Takes the natural logarithm of the input.
is.function(log) ## Checks whether the object is a function.
log(exp(1)) ## Note the exponential function written as "exp().
c(1, 2, 3, 4) ## Concatenate function which created a vector.

## ----Function.Parts.1----------------------------------------------------
str(paste)
formals(paste) ## Prints the arguments of a function.
body(paste) ## Prints the body of a function.
environment(paste) ## Prints the environment of a function.

getMethod("log") ## Shows the source code of one of the functions with the same name in the global environment. If this function does not work, try the below functions.
str(log) ## str function gives the structure of the function with its arguments. I rarely use this function for functions but it might be usefull to reveal the full structure of a function.

getAnywhere(log) ## Shows the information of matching function names in all pacakges.
getAnywhere(log)[1] ## Selecting the first function.

getAnywhere(paste) ## Since there is only one function with the mathing name, functions source code is revealed immediately.

getAnywhere(Head.Tail) ## This is a user-written R function. Note the curly braces which represents the body of the function.
formals(Head.Tail)
body(Head.Tail)
environment(Head.Tail)

# edit(log) ## Use "edit" function to open the source code of a function in a small editor window in RStudio.

## ----Function.Arguments.1------------------------------------------------
args(c) ## No argument exists.
args(log) ## Has one argument.
args(setdiff) ## Has two arguments.
args("+") ## Has two arguments.
args(mean) ## At least one argument. "..." means that some other arguments can be passed to other functions.
formals(mean)
args(nb2mat) ## Has 4 arguments and the last three have default values.
formals(nb2mat)

## ----Argument.Matching.1-------------------------------------------------
a <- c(1, 2, 3, 4) ## Vector 1.
b <- c(1, 2, 5, 9) ## Vector 2.

# setdiff: Everything in "x" and not in "y".
getAnywhere(setdiff)[3]
str(setdiff)

setdiff(x = a, y = b) ## With arguments.
setdiff(y = b, x = a) ## With arguments.
setdiff(a, b) ## Without arguments.
setdiff(b, x = a) ## With some arguments.

setdiff(b, a)

## ----Creating.R.Functions.Syntax.1, error = TRUE-------------------------
# R function syntax.
function.syntax <- function(arguments) {
    ## Do something interesting
}
str(function.syntax)

function.syntax(arguments) ## Calling the function.

## ----My.R.Functions.1----------------------------------------------------
# Function with no arguemnts.
myfunction <- function() {
    x <- rnorm(100)
    mean(x) ## The last expression will be returned.
}
myfunction()
formals(myfunction)
body(myfunction)
environment(myfunction)

# Function with one argument.
my.cube.func <- function(x) {
    x^3
}
my.cube.func(2)
formals(my.cube.func)
body(my.cube.func)
environment(my.cube.func)

# Function with one argument and return function.
my.variance <- function(x) { ## Sample variance.
    a <- (sum(x^2) - length(x) * mean(x)^2) / (length(x) - 1)
    return(a) ## You can also use return function to return a specific value.
}
my.variance(c(1:5))
var(c(1:5)) ## Built-in R function gives the same result.

## ----My.R.Functions.2, eval = FALSE--------------------------------------
## # R code chunk is not evaluated.
## 
## # Naming a user created R function with a known built-in R function name.
## c <- function(x) {
##     a <- x^x
##     return(a)
## }
## c(4)
## class(c)
## rm(c) ## The created c function is removed.

## ----Default.Values.1----------------------------------------------------
# Function with default values.
sum.of.squares <- function(x, About = mean(x)) { ## Sum of squares.
    x <- x[!is.na(x)]
    a <- sum((x - About)^2)
    a <- return(a)
}
sum.of.squares(x = c(-2:2))
sum.of.squares(x = c(-2:2), About = 0)
sum.of.squares(x = c(-2:2), About = 1)

# Function with default values including NULL value.
my.function <- function(x, Power = 2, Addition = 10, Remove.NA = NULL) { ## Just a function.
    a <- sum(x + Addition, na.rm = Remove.NA)^Power
    a <- return(a)
}
my.function(x = c(1:5), Power = 3, Addition = 1, Remove.NA = TRUE)
my.function(x = c(1:5, NA), Power = 3, Addition = 1, Remove.NA = FALSE)
my.function(x = c(1:5, NA), Power = 3, Addition = 1, Remove.NA = NULL)

## ----Nested.Functions.1--------------------------------------------------
# First function.
my.square <- function(x) {
    x <- x[!is.na(x)] ## NA values are subsetted.
    a <- x^2
    return(a)
}

# Second function.
my.cube <- function(x) {
    x <- x[!is.na(x)] ## NA values are subsetted.
    a <- x^3
    return(a)
}

# The main function which nests the first and second functions.
sum.square.cube <- function(x) {
    a <- sum(my.square(x))
    b <- sum(my.cube(x))
    c <- a + b
    return(c)
}
sum.square.cube(1)
sum.square.cube(2)
sum.square.cube(c(1:5))
sum.square.cube(c(1:5, rep(NA, 3)))

## ----Nested.Functions.2--------------------------------------------------
# Function which creates functions as output.
make.power <- function(power) {
    power.func <- function(base) {
        return(base^power)
        }
    return(power.func)
}
make.power(2)

# New function 1.
square.func <- make.power(2) ## Created a square function.
square.func(3) ## Takes the square of 2.

## New function 2.
cube.func <- make.power(3) ## Created a cube function.
cube.func(3) ## Takes the cube of 3.

# What's in a function's environment?
ls(environment(cube.func)) ## Gives the defined object in cube.func.
get("power", environment(cube.func)) ## Gives the value of "power" in cube.func.

## ----Scoping.Case.Study.1.1----------------------------------------------
# Case Study 1
y <- 10 ## A variable defined in global environment not inside of a function.

func.1 <- function(x) { ## Function 1.
    x*y
}

func.2 <- function(x) { ## Function 2.
    y <- 2 ## A variable defined in the environment of func.2 function.
    y^2 + func.1(x)
}

## Which y values does func.1 and func.2 use?
func.2(3) ## Check the result.
### With lexical scoping the value of y in the function func.1 is looked up in the environment in which the function is created, in this case the global environment, so the value of y is 10.
### For func.2() function, y value is 2 which is defined while the func.2 is created.

## ----Scoping.Case.Study.2.1----------------------------------------------
# Case Study 2
y <- 10 ## A variable defined in global environment not inside of a function.

func.2 <- function(x) { ## Function 2.
    y <- 2 ## A variable defined in the environment of func.2 function.

    func.1 <- function(x) { ## Function 1.
        x*y
    }
    
    y^2 + func.1(x)
}

## Which y values does func.1 and func.2 use?
func.2(3) ## Check the result.
### This time since func.1 is created inside the func.2, the func.1 will use 2 as the y value.
### Func.2 also uses value 2 for y.

## ----Conditional.Statements.Syntax.1, eval = FALSE-----------------------
## # R code chunk is not evaluated.
## 
## # Syntax Case 1
## if(conditional.statement) {
##     ## Executes the code if the conditional.statement is TRUE.
## }
## 
## # Syntax 2
## if(conditional.statement) {
##     ## Executes the code if the conditional.statement is TRUE.
## } else {
##     ## Executes the code if the conditional.statement is FALSE.
## }
## 
## # Syntax 3
## if(conditional.statement.1) {
##     ## Executes the code if the "conditional.statement.1" is TRUE.
## } else if (conditional.statement.2) {
##     ## Executes the code if the "conditional.statement.1" is FALSE but "conditional.statement.2" is TRUE.
## } else {
##     ## Executes the code if the "conditional.statement.1" and "conditional.statement.2" are FALSE.
## }
## 
## # Syntax 4
## if(conditional.statement.1) {
##     ## Executes the code if the "conditional.statement.1" is TRUE.
##     if (conditional.statement.2) {
##         ## Executes the code if the "conditional.statement.1" and "conditional.statement.2" is TRUE.
##             if (conditional.statement.2) {
##             ## Executes the code if the "conditional.statement.1", "conditional.statement.2" and "conditional.statement.2" are TRUE.
##             }
##     }
## }

## ----Conditional.Statements.Case.1, error = TRUE-------------------------
# Simple if (single) statement.
x <- -8
if (x < 0) {
   print("Input is a negative number.")
}

## Simple if (multiple) statements with message function.
x <- sample(x = c(-1000:1000), size = 1, replace = TRUE, prob = NULL)
if (x < 0) {
    print("Input is a negative number.")
    message(paste0("Your input is ", x)) ## You can use "message" function to give informational note on the console.
    warning("Something unusual is going on.")
}
if (x > 0) {
    print("Input is a positive number.")
    message(paste0("Your input is ", x))
    warning("Something unusual is going on.")
}

# Simple if (single) statement with stop function.
if (!("$" %in% letters)) { ## See how to use "stop" function to end the conditional statement if
    stop("Invalid letter.") ## Invalid letter.
}

## ----Conditional.Statements.Case.2---------------------------------------
# Simple if and else statemenst.
x <- sample(x = c(-1000:1000), size = 1, replace = TRUE, prob = NULL)
if (x < 0) {
   print(paste0("Input, ", x , ", is a negative number."))
} else {
   print(paste0("Input, ", x , ", is a negative number."))
}

# Simple if and else statement with value assigning inside the conditional statement.
a <- 10
if (a > 3) {
    b <- 10 ## Creating a new variable.
} else {
    b <- 0
}
b

# Simple if and else statement with value assigning.
x <- 50
y <- if (x == 50) {
    0
    } else {
        1
    }
y

## ----Conditional.Statements.Case.3---------------------------------------
# If, else if (single) and else statements.
x <- 5
if (x < 0) {
    print("x is a negative number.")
} else if (x == 0) {
    print("x is zero.")
} else {
    print("x is a positive number.")
}

# If, else if (multiple) and else statements.
grade <- 100
if (grade < 70) {
    print("Keep studying!!!")
} else if (grade < 80) {
    print("Average")
} else if (grade < 90) {
    print("Good")
} else if (grade < 100) {
    print("Very Good")
} else {
    print("Excellent")
}

## ----Conditional.Statements.Case.4---------------------------------------
# If-else statements are nested in if-else statements.
## This conditional statement yields the same answer.
grade <- 100
if (grade < 100) {
    if (grade < 90) {
        if (grade < 80) {
            if (grade < 70) {
                print("Keep studying!!!")
            } else {
                print("Average")
            }
        } else {
            print("Good")
        }
    } else {
        print("Very Good")
    }
} else {
    print("Excellent")
}

# If statements are nested in if statements.
## This conditional statement yields the same answer.
grade <- 100
if (grade < 100) {
    if (grade < 90) {
        if (grade < 80) {
            if (grade < 70) {
                print("Keep studying!!!")
            } 
            if (grade >= 70) {
                print("Average")
            }
        } 
        if (grade >= 80) {
            print("Good")
        }
    } 
    if (grade >= 90) {
        print("Very Good")
    }
} else {
    print("Excellent")
}

## ----ifelse.Function.1---------------------------------------------------
# If else statement for vectorized objects.
x <- 1:10
if (x < 5) {
    x <- 0
}
x

# ifelse function for vectorized object.
x <- 1:10
y <- ifelse(x < 5, 0, 1) 
y

## ----for.Loops.1, eval = FALSE-------------------------------------------
## # R code chunk is not evaluated.
## 
## # The general syntax of a for loop
## for (counter in sequence) {
##     ## Executes the code for each iteration of counter.
## }
## 
## # The most common syntax of a for loop
## for (i in x) {
##     ## Executes the code for each iteration of counter i in x.
## }

## ----for.Loops.2---------------------------------------------------------
# Simple for loop.
## This loop takes the variable "i" and for each iteration of the loop 1, 2, 3, ..., 10, are assigned to it. After the last iteration the loop exits.
for (i in 1:10) { ## Counter is "i".
    print(i)
}

# Simple for loop with a different counter.
x <- c("a", "b", "c")
for (NCSU in 1:length(x)) { ## Counter is "NCSU"
    print(x[NCSU])
}
NCSU

# Simple for loop with seq-along function.
sample.size <- sample(x = c(5:10), size = 1)
x <- sample(x = c(-10:10), size = sample.size, replace = TRUE, prob = NULL)
for (j in seq_along(x)) { ## Counter is "j"
    print(x[j])
}
j

## ----for.Loops.3---------------------------------------------------------
# Simple nested for loops.
x <- matrix(data = c(1:6), nrow = 2, ncol = 3)
x
for (i in 1:nrow(x)) { ## Looping over rows.
    for (j in 1:ncol(x)) { ## Looping over columns.
        print(x[i, j])
    }
}
i ## Number of rows.
j ## Number of columns.

## ----for.Loops.4---------------------------------------------------------
# Next statement.
for (i in 1:7) {
    if (i <= 5) {
        next ## Skip the first 5 iterations.
    }
    print(i)
}

# Break statement.
for (i in 1:7) {
    if (i > 5) {
        break ## Terminates the loop on the 6th iteration.
    }
    print(i)
}

## ----while.Loops.1-------------------------------------------------------
# Simple while loop
count <- 0 ## "Count" variable initializes with 0.
while (count < 10) {
    print(count)
    count <- count + 1 ## Count variable is updated and the loop starts again.
}

# Sometimes there will be more than one condition in the test.
x <- 5
while (x >= 3 & x <= 10) {
    print(x)
    coin <- rbinom(n = 1, size = 1, prob = 0.5) ## Flips a fair coin. 0 means fail, 1 means success.
    if (coin == 1) { ## random walk
        x <- x + 1
    } else {
        x <- x - 1
    }
}

## ----repeat.Loops.1------------------------------------------------------
# Simple repeat loop.
x <- 1 ## Initial value.
repeat {
   print(x)
   if (x == 6) {
       break ## If the condition TRUE then stop.
   } else {
       x <- x + 1 ## If the condition is FALSE then run the this code.
   }
}

## ----repeat.Loops.2, eval = FALSE----------------------------------------
## # R code chunk is not evaluated.
## 
## # Simple repeat loop which does not stop.
## x <- 1 ## Initial value.
## repeat {
##    print(x)
##    if (x > Inf) {
##        break ## If the condition TRUE then stop.
##    } else {
##        x <- x + 1 ## If the condition is FALSE then run the this code.
##    }
## }

## ----lapply.1------------------------------------------------------------
# lapply function with a list.
x <- list(a = 1:4, b = rnorm(10), c = rnorm(20, 1), d = rnorm(100, 5))
lapply(X = x, FUN = mean) ## Gives the mean of each element on a list.

# lapply function with a numeric vector
x <- c(1:4)
lapply(x, runif) ## "runif" function creates uniform rondom variables. First arguement in "runif" is the number of the variables that you want to create uniform random variables. lapply gives you runif(1), runif(2).... Note that "runif" has other arguments but we dont need to specify these right now since they have default values. The default is uniform between 0 and 1.

# lapply function with a numeric vector and passing arguments from other functions.
x <- c(1:4)
lapply(x, runif, min = 0, max = 10) ## "min" and "max" arguements are passed from "runif" function.

## ----lapply.2------------------------------------------------------------
# lapply function with a anonymous function.
x <- list(a = matrix(1:8, 4, 2), b = matrix(1:12, 3, 4))
lapply(x, function(col) col[, 1]) ## An anonymous function for extracting the first column of each matrix. There is no function "col" but we just write it and used in lapply. After lapply is finished this function will go away so this "elt" function is anonymous function.

## ----lapply.3------------------------------------------------------------
# Using split and lapply functions together.
x <- c(rnorm(5), runif(5), rnorm(5, 1))
f <- gl(n = 3, k = 5) ## Generates factor levels.
split(x, f)
lapply(split(x, f), mean)

## ----sapply.1------------------------------------------------------------
# sapply function with a list.
x <- list(a = 1:4, b = rnorm(10), c = rnorm(20, 1), d = rnorm(100, 5))
lapply(X = x, FUN = mean) ## List format.
sapply(X = x, FUN = mean, simplify = FALSE) ## Same as lapply.
sapply(X = x, FUN = mean) ## Vector format.
mean(x) ## Note that mean function cannot handle list objects.


## ----apply.1-------------------------------------------------------------
# apply function on a matrix which returns a vector.
x <- matrix(rnorm(200), 20, 10)
y <- apply(X = x, MARGIN = 2, FUN = mean) ## Means of columns.
y
class(y)
str(y)

apply(x, 1, sum) ## Calculates the sum of each row.

# apply function on a matrix which returns a array.
x <- matrix(rnorm(200), 20, 10)
y <- apply(X = x, MARGIN = 2, FUN = quantile, probs = c(0.25, 0.75)) ## Gives the first and the third quantiles of each column.
y
class(y)
str(y)

# apply function on an array which returns a matrix.
## Gives averages of an array in a matrix format.
x <- array(rnorm(2 * 2 * 10), c(2, 2, 10)) ## This array has 3 dimensions: with 2 rows, 2 columns and the 3rd dimension with number 10.
apply(x, c(1, 2), mean) ## Generates the mean of the array with the 1st and 2nd dimension. In other meaning, 3rd dimension is collapssed. So the resulting matrix will be a 2x2 matrix with means.
rowMeans(x, dims = 2) ## Gives the same result as above. "2" represents the first number of dimensions which are preserved.

# apply function on an array which returns an array.
x <- array(rnorm(2 * 2 * 10), c(2, 2, 10)) 
apply(x, c(2, 3), mean) ## Means for the 2nd and the 3rd dimensions.

## ----apply.2-------------------------------------------------------------
x <- matrix(rnorm(15), 3, 5)
apply(x, 1, sum) ## Same as "rowSums(x)" function.
rowSums(x)

apply(x, 1, mean) ## Same as "rowMeans(x)" function.
rowMeans(x)

apply(x, 2, sum) ## Same as "colSums(x)" function.
colSums(x)

apply(x, 2, mean) ## Same as "colMeans(x)" function.
colMeans(x)

## ----tapply.1------------------------------------------------------------
# tapply function (simple).
x <- c(rnorm(10), runif(10), rnorm(10, 1))
f <- gl(n = 3, k = 10) ## Generates factor levels.
tapply(X = x, INDEX = f, FUN = mean) ## A factor level is assigned to each value in x in order.
tapply(x, f, range) ## Gives the min and max within the subset of x.

tapply(x, f, mean, simplify = FALSE) ## The result is in a list.
lapply(split(x, f), mean) ## Same as above.

# tapply function (complex).
x <- c(rnorm(5), rnorm(5, 1), rnorm(5, 2), rnorm(5, 3)) ## Our values.
f1 <- factor(rep(1:2, each = 10)) ## First factor.
f2 <- factor(rep(rep(3:4, each = 5), times = 2)) ## Second factor.
f <- list(f1, f2) ## List of factors.
tapply(X = x, INDEX = f, FUN = mean)

## ----mapply.1------------------------------------------------------------
# mapply function (simple).
list(rep(1, 4), rep(2, 3), rep(3, 2), rep(4, 1)) # Instead we can do the below code.
mapply(rep, 1:4, 4:1) ## mapply function takes the arguement in order.

# mapply function (complex).
noise <- function(n, mean, sd) { ## A function for n, mean and sd.
    rnorm(n, mean, sd)
}
noise(5, 1, 2)
noise(1:5, 1:5, 2) ## It does not work correctly for set of n's and means's. No vectorization.
mapply(noise, 1:5, 1:5, 2) ## With mapply, it will be vectorized.

## ----Case.Study.Conditional.Statements.1---------------------------------
n <- 10 ## Number of employees.
all.incomes <- sample(x = c(0:(5*10^5)), size = n, replace = FALSE, prob = NULL) ## Random sample for income.
all.incomes ## Income values.

## ----Case.Study.Conditional.Statement.1----------------------------------
# Calculating for one employee.
income <- all.incomes[1] ## Income for the first employee.
if (income <= 0) {
    tax <- 0
} else if (income <= 9325) {
    tax <- income * 0.1
} else if (income <= 37950) {
    tax <- income * 0.15
} else if (income <= 91900) {
    tax <- income * 0.25
} else if (income <= 191650) {
    tax <- income * 0.28
} else if (income <= 416700) {
    tax <- income * 0.33
} else if (income <= 418400) {
    tax <- income * 0.35
} else {
    tax <- income * 0.396
}
c("Income" = income, "Tax" = tax, "Net Income" = income - tax)

## ----Case.Study.Loop.1---------------------------------------------------
# Calculating for one employee.
tax.brackets <- list(c(0, 0.1), c(9325, 0.15), c(37950, 0.25), c(91900, 0.28), c(191650, 0.33), c(416700, 0.35), c(418400, 0.396)) ## Tax brackets in a list.
income <- all.incomes[1]
for (i in 1:length(tax.brackets)) {
    if (tax.brackets[[i]][1] < income) {
        tax.rate <- tax.brackets[[i]][2]
        tax <- income * tax.rate
    }
}
c("Income" = income, "Tax Rate" = tax.rate, "Tax" = tax, "Net Income" = income - tax)

## ----Case.Study.Loop.2---------------------------------------------------
# Calculating for all employee.
tax.brackets <- list(c(0, 0.1), c(9325, 0.15), c(37950, 0.25), c(91900, 0.28), c(191650, 0.33), c(416700, 0.35), c(418400, 0.396)) ## Tax brackets in a list.

for (j in 1:length(all.incomes)) {
    income <- all.incomes[j]
    for (i in 1:length(tax.brackets)) {
        if (tax.brackets[[i]][1] < income) {
            tax.rate <- tax.brackets[[i]][2]
            tax <- income * tax.rate
        }
    }
    if (j == 1) {
        results <- c(income, tax.rate, tax, income - tax)
    } else {
        temp <- c(income, tax.rate, tax, income - tax)
        results <- rbind(results, temp)
    }
}
results <- as.data.frame(results, row.names = paste("Employee", " ", 1:length(all.incomes)), stringsAsFactors = FALSE)
colnames(results) <- c("Income", "Tax Rate", "Tax", "Net Income")
results

## ----Case.Study.Function.1-----------------------------------------------
# Calculating with any pre-specified tax brackets for all employees.
## Pre-specified tax brackets in a list.
tax.brackets <- list(c(0, 0.1), c(9325, 0.15), c(37950, 0.25), c(91900, 0.28), c(191650, 0.33), c(416700, 0.35), c(418400, 0.396))

## Function with Income and Tax.Brakets options.
tax.func <- function(Income, Tax.Brackets) {
    for (j in 1:length(Income)) {
        for (i in 1:length(Tax.Brackets)) {
            if (Tax.Brackets[[i]][1] < Income[j]) {
                tax.rate <- Tax.Brackets[[i]][2]
                tax <- Income[j] * tax.rate
            }
        }
        if (j == 1) {
            results <- c(Income[j], tax.rate, tax, Income[j] - tax)
        } else {
            temp <- c(Income[j], tax.rate, tax, Income[j] - tax)
            results <- rbind(results, temp)
        }
    }
    results <- as.data.frame(results, row.names = paste("Employee", " ", 1:length(Income)), stringsAsFactors = FALSE)
    colnames(results) <- c("Income", "Tax Rate", "Tax", "Net Income")
    return(results)
}

tax.func(Income = all.incomes, Tax.Brackets = tax.brackets)

## ----Case.Study.Function.2-----------------------------------------------
# Calculating with any pre-specified tax brackets for all employees.
## Pre-specified tax brackets in a list.
tax.brackets <- list(c(0, 0.1), c(50000, 0.25), c(100000, 0.30), c(150000, 0.35), c(200000, 0.40), c(300000, 0.45), c(400000, 0.50))

tax.func(Income = all.incomes, Tax.Brackets = tax.brackets)

## ----Session.Info, collapse = FALSE--------------------------------------
R.version.string ## Returns the R version in a string.
sessionInfo() ## From utils package.
devtools::session_info() ## From devtools package.

## ----Used.R.Functions.1, echo = FALSE------------------------------------
if (file.exists(paste0(objects.path, "used.r.functions.RData"))) {
    used.r.functions <- readRDS(file = paste0(objects.path, "used.r.functions.RData"), refhook = NULL) ## Loads the used.r.functions.RData file as a new object.
    used.r.functions
}

