proj.dir = file.path('~/masterarbeit_code')
fit.dir = file.path( proj.dir, 'fit' )
models.dir = file.path( proj.dir, 'models' )
dataprep.dir = file.path( proj.dir, 'dataprep')
data.dir = file.path( dataprep.dir, 'dat')
post.dir = file.path( proj.dir, 'posterior' )
separator = .Platform$file.sep

assertthat::assert_that(dir.exists(proj.dir),
                        dir.exists(models.dir),
                        dir.exists(dataprep.dir),
                        dir.exists(data.dir),
                        dir.exists(post.dir),
                        dir.exists(fit.dir)
                        )
