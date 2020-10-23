$ontext


    @author: Torbjörn Jansson (SLU)

$offtext

*===============================================================================
* Initialization
*===============================================================================


* --- Create useful environment variables

$setGlobal RESULTS_OUT "."
$setLocal ERROR_FILE "%RESULTS_OUT%\error_3_countries.gdx"


* --- Declare symbols to use

set s_reg "List of regions, not yet a set" /eu,usa,china/;
set regtot "Regions and totals" /#s_reg, na "Not applicable here", Total "Total"/
set reg(regtot) "Regions in the model" / #s_reg/;

alias(reg,reg1,reg2);

set dataType "Types of data in the model" /raw, cal, sim/;

set cols "Columns in database" /
    pConsDem "Average consumer price of domestic goods and imports"
    pConsAgg "Consumer price index of aggregate good (price of utility)"
    pConsDom "Consumer price of domestic good"
    pConsImp "Consumer price of imported good"
    pImpAgg  "Price index of imported good (price of utility)"

    pprod     "Producer price"
    pfob      "Free-On-Board price"
    pcif      "Price with Cost, Insurance, and Freigh in importing country"
    ppmrg     "Producer price markup for domestic market"


    qSup      "Quantity produced"
    qDom      "Quantity demanded that is of domestic origin"
    qImp      "Quantity demanded that is imported, sum of all flows"
    qDem      "Quantity demanded (domestic plus imports)"
    uDem      "Utility of aggregated demand on top level"
    uImp      "Utility of import aggregate"
    qExp      "Quantity exported from region"
    qTrd      "Quantity of import flow, physical units, bilateral data"

    tc        "Unit trade cost"

    /;

set rows "Rows in database" /beef/;
set comm(rows) "Commodities" /beef/;

set problemData(regtot,regtot,cols,rows,datatype) "Data problems detected";
parameter p_problemDiff(regtot,regtot,cols,rows,datatype) "Size of data problem, diff";
scalar p_problemTol "Tolerance for absolute deviation" /1E-6/;
set colsToCheck(cols) "Columns to check something for";


table p_marketData(regtot,cols,rows,dataType) Data for each market (no bilateral data)
                            raw
eu.qsup.beef                80
usa.qsup.beef               125
china.qsup.beef             85

eu.pprod.beef               1
usa.pprod.beef              1
china.pprod.beef            1

eu.pConsDom.beef            2
usa.pConsDom.beef           2
china.pConsDom.beef         2

eu.qdom.beef                70
usa.qdom.beef               85
china.qdom.beef             80

eu.qimp.beef                20
usa.qimp.beef               15
china.qimp.beef             20
;


table p_tradeData(regtot,regtot,cols,rows,dataType) "Bilateral data relating to trade"
                    raw
eu.usa.qTrd.beef     20
eu.china.qTrd.beef    0
usa.eu.qTrd.beef     10
usa.china.qTrd.beef   5
china.eu.qTrd.beef    0
china.usa.qTrd.beef  20
;


parameter p_tc(reg,reg1) "Unit trade costs for import of reg from reg1";
p_tc(reg,reg1) = 1;

parameter p_report(regtot,regtot,cols,rows,datatype);




*-------------------------------------------------------------------------------
* Declare the economic model
*-------------------------------------------------------------------------------

parameter p_scaleTop(reg,comm) "CES scale parameter for top nest";
parameter p_alphaTop(reg,comm) "CES share parameter for top nest";
parameter p_subsElasTop(reg,comm) "Elasticity of substitution in top nest";

parameter p_scaleImp(reg,comm) "CES scale parameter for import nest";
parameter p_alphaImp(reg,reg,comm) "CES share parameter for import nest";
parameter p_subsElasImp(reg,comm) "Substitution parameter in import nest";

parameter p_supScale(reg,comm) "Supply function scale parameter";
parameter p_supElas(reg,comm) "Supply elasticity";


* --- Variables in the model

variable v_qDom(reg,comm) "Consumption of domestically produced good";
variable v_qImp(reg,comm) "Consumption of imported good";
variable v_qTrd(reg,reg1,comm) "Trade defined as import of reg from reg1";
variable v_qSup(reg,comm) "Supply quantity";
variable v_armImpPriceIndex(reg,comm) "Price index of import utility aggregate";

variable v_pConsDom(reg,comm) "Consumer price of domestically produced good";

positive variables v_qDom,v_qImp, v_pConsDom;


* --- Equations in the model

equation e_armSharesTop(reg,comm) "Allocation of demand to domestic and import sources";
equation e_armAggTop(reg,comm) "Quantity restriction on aggregated goods";
equation e_armImpDemand(reg,reg,comm) "Allocation of total import demand to each origin";
equation e_armImpPriceIndex(reg,comm) "Price index associated with the utility of the import aggregate";

equation e_qSup(reg,comm) "Supply function";
equation e_pConsDom(reg,comm) "Relating the consumer price of domestic good to domestic producer price";


* --- Define the ratio of domestic to imports

e_armSharesTop(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    v_qDom(reg,comm) / v_qImp(reg,comm)
        =E=
    [p_marketData(reg,"pConsDom",comm,"cal")/p_marketData(reg,"pConsImp",comm,"cal")
    *(1-p_alphaTop(reg,comm))/(p_alphaTop(reg,comm))]**(-p_subsElasTop(reg,comm));


* --- Closing the balance at the top: constrain the utility level (separability!)

e_armAggTop(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    p_marketData(reg,"uDem",comm,"cal")
        =E=

    p_scaleTop(reg,comm)
    * [
          p_alphaTop(reg,comm)
        * v_qDom(reg,comm)**((p_subsElasTop(reg,comm)-1)/p_subsElasTop(reg,comm))

        + (1-p_alphaTop(reg,comm))
        * v_qImp(reg,comm)**((p_subsElasTop(reg,comm)-1)/p_subsElasTop(reg,comm))

      ]**(p_subsElasTop(reg,comm)/(p_subsElasTop(reg,comm)-1));


* --- Second level nest: define the quantity demanded from different sources
*     based on minimization of cost of reaching utility uImp

e_armImpDemand(reg,reg1,comm) $ p_tradeData(reg,reg1,"qTrd",comm,"cal")..
    v_qTrd(reg,reg1,comm)
        =E=
    p_alphaImp(reg,reg1,comm)**p_subsElasImp(reg,comm)
    * p_scaleImp(reg,comm)**(p_subsElasImp(reg,comm)-1)
    * (v_armImpPriceIndex(reg,comm)/p_tradeData(reg,reg1,"pCif",comm,"cal"))**p_subsElasImp(reg,comm)
    * p_marketData(reg,"uImp",comm,"cal");

e_armImpPriceIndex(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    v_armImpPriceIndex(reg,comm)
        =E=
    1/p_scaleImp(reg,comm)
    *sum[reg1 $ p_tradeData(reg,reg1,"qTrd",comm,"cal"),
            p_alphaImp(reg,reg1,comm)**p_subsElasImp(reg,comm)
          * p_tradeData(reg,reg1,"pCif",comm,"cal")**(1-p_subsElasImp(reg,comm))
        ]**(1/(1-p_subsElasImp(reg,comm)));


* --- Put supply in relation to prices (a constant elasticity function without cross effects)

e_qSup(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    v_qSup(reg,comm)
        =E= p_supScale(reg,comm)*p_marketData(reg,"pprod",comm,"cal")**p_supElas(reg,comm);


* --- Price linkage equations

e_pConsDom(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    v_pConsDom(reg,comm) =E= p_marketData(reg,"pProd",comm,"cal") + p_marketData(reg,"ppMrg",comm,"cal");



model m_trade "Model of trade" /e_armSharesTop, e_armAggTop
                                e_qSup
                                e_armImpDemand, e_armImpPriceIndex
                                e_pConsDom/;




*===============================================================================
* Data preparation and calibration
*===============================================================================

*-------------------------------------------------------------------------------
* Preparation of data
* Balance data in calibration point based on RAW data
*-------------------------------------------------------------------------------

p_marketData(reg,cols,comm,"cal") = p_marketData(reg,cols,comm,"raw");

p_marketData(reg,"qdem",comm,"cal")
    = p_marketData(reg,"qdom",comm,"cal") + p_marketData(reg,"qimp",comm,"cal");

p_marketData(reg,"qdem",comm,"cal")
    = p_marketData(reg,"qdom",comm,"cal") + p_marketData(reg,"qimp",comm,"cal");

p_marketData(reg,"qExp",comm,"cal")
    = p_marketData(reg,"qSup",comm,"cal") - p_marketData(reg,"qDom",comm,"cal");


* --- Assume for now that import prices are equal domestic prices

p_marketData(reg,"pConsImp",comm,"cal") = p_marketData(reg,"pConsDom",comm,"cal");


* --- The average consumer price exhausts the budget

p_marketData(reg,"pConsDem",comm,"cal")
    = (  p_marketData(reg,"pConsDom",comm,"cal")*p_marketData(reg,"qdom",comm,"cal")
       + p_marketData(reg,"pConsImp",comm,"cal")*p_marketData(reg,"qimp",comm,"cal"))
      / p_marketData(reg,"qimp",comm,"cal");


* --- The domestic consumer price is the domestic producer price plus mark-up

p_marketData(reg,"ppmrg",comm,"cal")
    = p_marketData(reg,"pConsDom",comm,"cal") - p_marketData(reg,"pprod",comm,"cal");


* --- Use RAW trade data to create balanced trade flows (now: assume it is balanced)

p_tradeData(reg,reg1,cols,rows,"cal") = p_tradeData(reg,reg1,cols,rows,"raw");


* --- Create FOB and CIF prices where there is trade

*   - Assume that the producer price margin is identical for exports as for domestic market?
p_tradeData(reg,reg1,"pFob",comm,"cal") $ p_tradeData(reg,reg1,"qTrd",comm,"cal")
    = p_marketData(reg1,"pProd",comm,"cal") + p_marketData(reg1,"ppmrg",comm,"cal");

p_tradeData(reg,reg1,"pCif",comm,"cal") $ p_tradeData(reg,reg1,"qTrd",comm,"cal")
    = p_tradeData(reg,reg1,"pFob",comm,"cal") + p_tc(reg,reg1);



*-------------------------------------------------------------------------------
* TESTS for data preparation part
* As many tests as possible to assert that the data is prepared as it should.
* "Only code if there is a failing test"
*-------------------------------------------------------------------------------

* --- Data preparation: if there is "raw" data, then there should be a "cal"
option kill=problemData;
option kill=colsToCheck;
colsToCheck("qdom") = yes;
colsToCheck("qimp") = yes;
colsToCheck("qdem") = yes;
colsToCheck("pConsDom") = yes;
colsToCheck("pConsImp") = yes;
colsToCheck("pConsDem") = yes;

problemData(reg,"na",colsToCheck,rows,"cal")
    $ [p_marketData(reg,colsToCheck,rows,"raw") and (not p_marketData(reg,colsToCheck,rows,"cal"))] = yes;

if(card(problemData),
    execute_unload "%ERROR_FILE%";
    abort "Raw demad data is there, but no CAL (calibrated) position", problemData;
);


* --- Data preparation: verify that certain accounting identities hold in calibration point
option kill=problemData;

problemData(reg,"na","qdem",comm,"cal")
    $ [p_marketData(reg,"qdem",comm,"cal")
            ne
        (p_marketData(reg,"qdom",comm,"cal") + p_marketData(reg,"qimp",comm,"cal"))] = yes;

if(card(problemData),
    execute_unload "%ERROR_FILE%";
    abort "Calibrated demand is not the sum of imports plus domestic", problemData;
);


* --- Assert that global production equals global consumption
option kill=problemData;

problemData("Total","na","qdem",comm,"cal")
    $ [sum(reg, p_marketData(reg,"qdem",comm,"cal"))
        ne
       sum(reg, p_marketData(reg,"qsup",comm,"cal"))] = yes;

if(card(problemData),
    execute_unload "%ERROR_FILE%";
    abort "Global supply does not equal global demand", problemData;
);


* --- Data preparation: verify that there are prices where there is a quantity
option kill=problemData;

*   - For total demand
problemData(reg,"na","pConsDem",comm,"cal")
    $ [p_marketData(reg,"qdem",comm,"cal") and (not p_marketData(reg,"pConsDem",comm,"cal"))] = yes;

*   - For domestically produced goods
problemData(reg,"na","pConsDom",comm,"cal")
    $ [p_marketData(reg,"qdom",comm,"cal") and (not p_marketData(reg,"pConsDom",comm,"cal"))] = yes;

*   - For imported goods
problemData(reg,"na","pConsImp",comm,"cal")
    $ [p_marketData(reg,"qimp",comm,"cal") and (not p_marketData(reg,"pConsImp",comm,"cal"))] = yes;

*   - That CIF price exists if there is trade
problemData(reg,reg1,"pCif",comm,"cal")
    $ [p_tradeData(reg,reg1,"qTrd",comm,"cal") and (not p_tradeData(reg,reg1,"pCif",comm,"cal"))] = yes;


$batinclude "assert_that_set_is_empty.gms" problemData "Price is missing albeit a quantity is there" %ERROR_FILE%


* --- The domestic consumer price is the domestic producer price plus mark-up
option kill=problemData;

problemData(reg,"na","pConsDem",comm,"cal")
    $ [p_marketData(reg,"ppmrg",comm,"cal")
         ne (p_marketData(reg,"pConsDom",comm,"cal")-p_marketData(reg,"pprod",comm,"cal"))] = yes;

$batinclude "assert_that_set_is_empty.gms" problemData "The domestic consumer price is not the producer price plus markup" %ERROR_FILE%


* --- Trade must balance for each region and commodity
option kill=problemData;

problemData(reg,"na","qImp",comm,"cal")
    $ [p_marketData(reg,"qImp",comm,"cal") ne sum(reg1, p_tradeData(reg,reg1,"qTrd",comm,"cal"))] = yes;

problemData(reg,"na","qExp",comm,"cal")
    $ [p_marketData(reg,"qExp",comm,"cal") ne sum(reg1, p_tradeData(reg1,reg,"qTrd",comm,"cal"))] = yes;


$batinclude "assert_that_set_is_empty.gms" problemData "Sum of trade flows do not match aggregate trade balance" %ERROR_FILE%


*-------------------------------------------------------------------------------
* Calibration of behavioural functions
*-------------------------------------------------------------------------------

* --- Calibrate top demand nest

p_marketData(reg,"uDem",comm,"cal") = p_marketData(reg,"qDem",comm,"cal");

p_marketData(reg,"pConsAgg",comm,"cal") = p_marketData(reg,"pConsDem",comm,"cal");

p_subsElasTop(reg,comm) $ p_marketData(reg,"qDem",comm,"cal") = 2;

p_alphaTop(reg,comm) $ p_marketData(reg,"qDem",comm,"cal")
    = [   p_marketData(reg,"pConsDom",comm,"cal")/p_marketData(reg,"pConsImp",comm,"cal")
       * (p_marketData(reg,"qDom",comm,"cal")/p_marketData(reg,"qImp",comm,"cal"))**(1/p_subsElasTop(reg,comm)) ]

    / [1+ p_marketData(reg,"pConsDom",comm,"cal")/p_marketData(reg,"pConsImp",comm,"cal")
       * (p_marketData(reg,"qDom",comm,"cal")/p_marketData(reg,"qImp",comm,"cal"))**(1/p_subsElasTop(reg,comm)) ]
    ;

p_scaleTop(reg,comm) $ p_marketData(reg,"qDem",comm,"cal")
    = p_marketData(reg,"qDem",comm,"cal")

    / [   p_alphaTop(reg,comm)
        * p_marketData(reg,"qDom",comm,"cal")**((p_subsElasTop(reg,comm)-1)/p_subsElasTop(reg,comm))

        + (1-p_alphaTop(reg,comm))
        * p_marketData(reg,"qImp",comm,"cal")**((p_subsElasTop(reg,comm)-1)/p_subsElasTop(reg,comm))

      ]**(p_subsElasTop(reg,comm)/(p_subsElasTop(reg,comm)-1));


* --- Calibrate supply functions

p_supElas(reg,comm) $ p_marketData(reg,"qDem",comm,"cal") = 1;

p_supScale(reg,comm)
    = p_marketData(reg,"qSup",comm,"cal")
    / p_marketData(reg,"pprod",comm,"cal")**p_supElas(reg,comm);


* --- Calibrate import CES functions

p_marketData(reg,"uImp",comm,"cal") $ p_marketData(reg,"qDem",comm,"cal")
    = p_marketData(reg,"qImp",comm,"cal");

p_marketData(reg,"pImpAgg",comm,"cal") $ p_marketData(reg,"qDem",comm,"cal")
    = p_marketData(reg,"pConsImp",comm,"cal");

p_subsElasImp(reg,comm) $ p_marketData(reg,"qDem",comm,"cal")  = 4;

parameter p_temp(reg,reg1,reg2,comm) "Intermediate computational result";

p_temp(reg,reg2,reg1,comm) $ [  p_tradeData(reg,reg2,"qTrd",comm,"cal")
                            and p_tradeData(reg,reg1,"qTrd",comm,"cal")
                            and (not sameas(reg2,reg1))]
    = p_tradeData(reg,reg2,"pCif",comm,"cal")
    / p_tradeData(reg,reg1,"pCif",comm,"cal")
    * [  p_tradeData(reg,reg2,"qTrd",comm,"cal")
       / p_tradeData(reg,reg1,"qTrd",comm,"cal")]**(1/p_subsElasImp(reg,comm));

p_alphaImp(reg,reg1,comm) $ p_tradeData(reg,reg1,"qTrd",comm,"cal")
    = 1 / [1 + sum(reg2 $ p_tradeData(reg,reg2,"qTrd",comm,"cal"), p_temp(reg,reg2,reg1,comm))];

p_scaleImp(reg,comm) $ p_marketData(reg,"qDem",comm,"cal")
    = p_marketData(reg,"uImp",comm,"cal")
    / sum[reg1, p_alphaImp(reg,reg1,comm)
               *p_tradeData(reg,reg1,"qTrd",comm,"cal")**((p_subsElasImp(reg,comm)-1)/p_subsElasImp(reg,comm))

         ]**(p_subsElasImp(reg,comm)/(p_subsElasImp(reg,comm)-1));


* --- Set starting values and bounds
v_qImp.L(reg,comm) $ p_marketData(reg,"qDem",comm,"cal") = p_marketData(reg,"qImp",comm,"cal");
v_qDom.L(reg,comm) $ p_marketData(reg,"qDem",comm,"cal") = p_marketData(reg,"qDom",comm,"cal");
v_qSup.L(reg,comm) $ p_marketData(reg,"qDem",comm,"cal") = p_marketData(reg,"qSup",comm,"cal");


*-------------------------------------------------------------------------------
* TESTS for successful calibration
* "Only code if there is a failing test"
*-------------------------------------------------------------------------------


* --- Calibration of demand: verify that the aggregate utility equals quantity
*     on the top level, for quantities and prices
option kill=problemData;

problemData(reg,"na","uDem",comm,"cal")
    $ [p_marketData(reg,"uDem",comm,"cal") ne p_marketData(reg,"qDem",comm,"cal")] = yes;

problemData(reg,"na","pConsAgg",comm,"cal")
    $ [p_marketData(reg,"pConsAgg",comm,"cal") ne p_marketData(reg,"pConsDem",comm,"cal")] = yes;

problemData(reg,"na","uImp",comm,"cal")
    $ [p_marketData(reg,"uImp",comm,"cal") ne p_marketData(reg,"qImp",comm,"cal")] = yes;

problemData(reg,"na","pImpAgg",comm,"cal")
    $ [p_marketData(reg,"pImpAgg",comm,"cal") ne p_marketData(reg,"pConsImp",comm,"cal")] = yes;

$batinclude "assert_that_set_is_empty.gms" problemData "Utility aggregate does not match quantity in calibration point" %ERROR_FILE%


* --- Check that Armington share parameters add up to unity
option kill=problemData;
option kill=p_problemDiff;

p_problemDiff(reg,"na","uImp",comm,"cal") = sum(reg1, p_alphaImp(reg,reg1,comm))-1;
problemData(reg,"na","uImp",comm,"cal") $ p_problemDiff(reg,"na","uImp",comm,"cal") = yes;

$batinclude "assert_that_set_is_empty.gms" problemData "Armington import share parameters do not add to unity." %ERROR_FILE%


* --- Just solving the model is a test in itself. Solve to verify calibration.

solve m_trade using cns;


* --- Calibration of demand: verify that simulated quantities
*     equal calibrated data counterparts
option kill=problemData;
option kill=p_problemDiff;

*   - For demand quantities
problemData(reg,"na","qDom",comm,"cal")
    $ [v_qDom.L(reg,comm) ne p_marketData(reg,"qDom",comm,"cal")] = yes;

*   - For supply quantities
problemData(reg,"na","qSup",comm,"cal")
    $ [v_qSup.L(reg,comm) ne p_marketData(reg,"qSup",comm,"cal")] = yes;

*   - For demand prices
problemData(reg,"na","pConsDom",comm,"cal")
    $ [v_pConsDom.L(reg,comm) ne p_marketData(reg,"pConsDom",comm,"cal")] = yes;

*   - For traded quantities
p_problemDiff(reg,reg1,"qTrd",comm,"cal")
    = v_qTrd.L(reg,reg1,comm) - p_tradeData(reg,reg1,"qTrd",comm,"cal");
problemData(reg,reg1,"qTrd",comm,"cal") $ [abs(p_problemDiff(reg,reg1,"qTrd",comm,"cal")) > p_problemTol] = yes;


$batinclude "assert_that_set_is_empty.gms" problemData "Some price or quantity is not calibrating correctly" %ERROR_FILE%



