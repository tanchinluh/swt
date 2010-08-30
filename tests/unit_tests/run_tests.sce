mode(-1);

//tmp = warning("query");
//warning("off");
//dir = getcwd();
chdir(get_absolute_file_path("run_tests.sce"));
seed = rand("seed");

printf("testing filter... ");
exec("filter_test.tst");
printf("ok!\n");

printf("testing utility functions... ");
exec("utility_test.tst");
printf("ok!\n");

printf("testing dwt1d ... ");
exec("dwt1d_test.tst");
printf("ok!\n");

printf("testing dwt2d ... ");
exec("dwt2d_test.tst");
printf("ok!\n");

printf("testing dwt2d 2.test... ");
exec("dwt2d_test2.tst");
printf("ok!\n");




rand("seed",seed);
//chdir(dir);
//warning(tmp);