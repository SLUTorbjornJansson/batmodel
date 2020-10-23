$ontext

    @author: Torbjörn Jansson, SLU
    @purpose: Verify that a set is empty, and abort with message if not

    @arguments: 1. "Problem Set" that should be empty
                2. String with error text to return
                3. full path and name of file to contain errors

    @date: 2020-10-22

$offtext

$setLocal PROBLEM_SET %1
$setLocal ERROR_MESSAGE %2
$setLocal ERROR_FILE %3

if(card(%PROBLEM_SET%),
    display "Test failure! All data unloaded to %ERROR_FILE%";
    execute_unload "%ERROR_FILE%";
    abort "%ERROR_MESSAGE%", %PROBLEM_SET%;
);