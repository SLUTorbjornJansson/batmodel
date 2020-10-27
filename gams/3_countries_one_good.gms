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
    pDemAgg  "Consumer price index of aggregate good (price of utility)"
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
variable v_qImp(reg,comm) "Consumption of imported good (physical quantity)";
variable v_uImp(reg,comm) "Consumer utility of import bundle";
variable v_qTrd(reg,reg1,comm) "Trade defined as import of reg from reg1";
variable v_qSup(reg,comm) "Supply quantity";
variable v_armImpPriceIndex(reg,comm) "Price index of import utility aggregate";

variable v_pProd(reg,comm) "Producer price";
variable v_pFob(reg,comm) "Export price (free-on-board) including consumer price margin";
variable v_pConsDom(reg,comm) "Consumer price of domestically produced good";

positive variables v_qDom,v_qImp, v_pConsDom;


* --- Equations in the model

equation e_armSharesTop(reg,comm) "Allocation of demand to domestic and import sources";
equation e_armAggTop(reg,comm) "Quantity restriction on aggregated goods";
equation e_armAggImp(reg,comm) "Aggregating import flows to obtain consumer utility of imports";
equation e_armImpDemand(reg,reg,comm) "Allocation of total import demand to each origin";
equation e_armImpPriceIndex(reg,comm) "Price index associated with the utility of the import aggregate";

equation e_supBalance(reg,comm) "Supply balance";
equation e_qSup(reg,comm) "Supply function";

equation e_pFob(reg,comm) "Export price definition";
equation e_pConsDom(reg,comm) "Relating the consumer price of domestic good to domestic producer price";


* --- Eq. 1: Define the ratio of domestic to imports

e_armSharesTop(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    v_qDom(reg,comm) / v_uImp(reg,comm)
        =E=
    [p_marketData(reg,"pConsDom",comm,"cal")/p_marketData(reg,"pConsImp",comm,"cal")
    *(1-p_alphaTop(reg,comm))/(p_alphaTop(reg,comm))]**(-p_subsElasTop(reg,comm));


* --- Eq. 2: Closing the balance at the top: constrain the utility level (separability!)

e_armAggTop(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    p_marketData(reg,"uDem",comm,"sim")
        =E=

    p_scaleTop(reg,comm)
    * [
          p_alphaTop(reg,comm)
        * v_qDom(reg,comm)**((p_subsElasTop(reg,comm)-1)/p_subsElasTop(reg,comm))

        + (1-p_alphaTop(reg,comm))
        * v_uImp(reg,comm)**((p_subsElasTop(reg,comm)-1)/p_subsElasTop(reg,comm))

      ]**(p_subsElasTop(reg,comm)/(p_subsElasTop(reg,comm)-1));


* --- Eq. 3: Second level nest: define the quantity demanded from different sources
*     based on minimization of cost of reaching utility uImp

e_armImpDemand(reg,reg1,comm) $ p_tradeData(reg,reg1,"qTrd",comm,"cal")..
    v_qTrd(reg,reg1,comm)
        =E=
    p_alphaImp(reg,reg1,comm)**p_subsElasImp(reg,comm)
    * p_scaleImp(reg,comm)**(p_subsElasImp(reg,comm)-1)
    * (v_armImpPriceIndex(reg,comm)/p_tradeData(reg,reg1,"pCif",comm,"cal"))**p_subsElasImp(reg,comm)
    * v_uImp(reg,comm);


* --- Eq. 4: Defining the Armington imports aggregate in terms of trade flows.
*            Not used, as the variable is implicitely defined by the top CES

e_armAggImp(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    v_uImp(reg,comm)
        =E=
    p_scaleImp(reg,comm)
    * SUM[reg1 $ p_tradeData(reg,reg1,"qTrd",comm,"cal"),
                p_alphaImp(reg,reg1,comm)
            *   v_qTrd(reg,reg1,comm)**((p_subsElasImp(reg,comm)-1)/p_subsElasImp(reg,comm))
          ]**(p_subsElasImp(reg,comm)/(p_subsElasImp(reg,comm)-1));


* --- Eq. 5: The Armington price index that is consistent with Utility for
*            weak separability.

e_armImpPriceIndex(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    v_armImpPriceIndex(reg,comm)
        =E=
    1/p_scaleImp(reg,comm)
    *sum[reg1 $ p_tradeData(reg,reg1,"qTrd",comm,"cal"),
            p_alphaImp(reg,reg1,comm)**p_subsElasImp(reg,comm)
          * p_tradeData(reg,reg1,"pCif",comm,"cal")**(1-p_subsElasImp(reg,comm))
        ]**(1/(1-p_subsElasImp(reg,comm)));


* --- Eq. 6: The market balance for supply: production equals domestic consumption
*            plus exports.

e_supBalance(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    v_qSup(reg,comm)
        =E=
    v_qDom(reg,comm)
    + sum(reg1 $ p_tradeData(reg1,reg,"qTrd",comm,"cal"), v_qTrd(reg1,reg,comm));


* --- Eq. 7: Put supply in relation to prices (a constant elasticity function without cross effects)

e_qSup(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    v_qSup(reg,comm)
        =E= p_supScale(reg,comm)*v_pProd(reg,comm)**p_supElas(reg,comm);


* --- Eq. 8: Linking consumer price of domestic good to domestic producer price

e_pConsDom(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    v_pConsDom(reg,comm) =E= v_pProd(reg,comm) + p_marketData(reg,"ppMrg",comm,"sim");


e_pFob(reg,comm) $ p_marketData(reg,"qdem",comm,"cal") ..
    v_pFob(reg,comm)
        =E= v_pProd(reg,comm) + p_marketData(reg,"ppMrg",comm,"sim");


model m_trade "Model of trade" /e_armSharesTop
*                                e_armAggImp
                                e_armAggTop
                                e_armImpDemand
                                e_armImpPriceIndex
                                e_supBalance
                                e_qSup
                                e_pConsDom
                                e_pFob/;




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


* --- The domestic consumer price is the domestic producer price plus mark-up

p_marketData(reg,"ppmrg",comm,"cal")
    = p_marketData(reg,"pConsDom",comm,"cal") - p_marketData(reg,"pprod",comm,"cal");


* --- Use RAW trade data to create balanced trade flows (now: assume it is balanced)

p_tradeData(reg,reg1,cols,rows,"cal") = p_tradeData(reg,reg1,cols,rows,"raw");


* --- Create FOB and CIF prices where there is trade
p_tradeData(reg,reg1,"pCif",comm,"cal") $ p_tradeData(reg,reg1,"qTrd",comm,"cal")
    = p_tradeData(reg,reg1,"pFob",comm,"cal")
    + p_tradeData(reg,reg1,"tc",comm,"cal");


*   - Assume that the producer price margin is identical for exports as for domestic market?
p_tradeData(reg,reg1,"pFob",comm,"cal") $ p_tradeData(reg,reg1,"qTrd",comm,"cal")
    = p_marketData(reg1,"pProd",comm,"cal") + p_marketData(reg1,"ppmrg",comm,"cal");


p_tradeData(reg,"pCif",comm,"cal")
    = p_marketData(reg,"pFob",comm,"cal")
    + p_tradeData(reg,"tc",comm,"cal");


p_marketData(reg,"pConsImp",comm,"cal") = p_marketData(reg,"pConsDom",comm,"cal");


* --- The average consumer price exhausts the budget

p_marketData(reg,"pConsDem",comm,"cal")
    = (  p_marketData(reg,"pConsDom",comm,"cal")*p_marketData(reg,"qdom",comm,"cal")
       + p_marketData(reg,"pConsImp",comm,"cal")*p_marketData(reg,"qimp",comm,"cal"))
      / p_marketData(reg,"qimp",comm,"cal");





*   - Define price index of imports to equal average import price in balanced point
p_marketData(reg,"pImpAgg",comm,"cal")
    = sum(reg1 $ p_tradeData(reg,reg1,"qTrd",comm,"cal"),
            p_tradeData(reg,reg1,"qTrd",comm,"cal")
           *p_tradeData(reg,reg1,"pCif",comm,"cal"))
    / p_marketData(reg,"qImp",comm,"cal");


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

problemData(reg,"na","pConsDom",comm,"cal")
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


* --- The price index of the import aggregate shall be the average import price
*
problemData(reg,"na","pImpAgg",comm,"cal")
    $ [p_marketData(reg,"pImpAgg",comm,"cal") ne p_marketData(reg,"pConsImp",comm,"cal")] = yes;

$batinclude "assert_that_set_is_empty.gms" problemData "Price index of imports does not match average import price." %ERROR_FILE%




*-------------------------------------------------------------------------------
* Calibration of behavioural functions
*-------------------------------------------------------------------------------

* --- Calibrate top demand nest

p_marketData(reg,"uDem",comm,"cal") = p_marketData(reg,"qDem",comm,"cal");

p_marketData(reg,"pDemAgg",comm,"cal") = p_marketData(reg,"pConsDem",comm,"cal");

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


* --- Compute the average import price. In this point (only!) it should equal
*     the import price index.

p_marketData(reg,"pImpAgg",comm,"cal")
    = sum(reg1 $ p_tradeData(reg,reg1,"qTrd",comm,"cal"),
            p_tradeData(reg,reg1,"qTrd",comm,"cal")
           *p_tradeData(reg,reg1,"pCif",comm,"cal"))
    / p_marketData(reg,"qImp",comm,"cal");


* --- Set starting values and bounds
v_uImp.L(reg,comm) $ p_marketData(reg,"qDem",comm,"cal") = p_marketData(reg,"uImp",comm,"cal");
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

problemData(reg,"na","pDemAgg",comm,"cal")
    $ [p_marketData(reg,"pDemAgg",comm,"cal") ne p_marketData(reg,"pConsDem",comm,"cal")] = yes;

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

p_marketData(reg,"uDem",comm,"sim") = p_marketData(reg,"uDem",comm,"cal")*1.0;
p_marketData(reg,"ppMrg",comm,"sim") = p_marketData(reg,"ppMrg",comm,"cal")*1.0;
solve m_trade using cns;


* --- Calibration of demand: verify that simulated quantities
*     equal calibrated data counterparts
option kill=problemData;
option kill=p_problemDiff;

*   - For demand quantities
p_problemDiff(reg,"na","qDom",comm,"cal")
    = v_qDom.L(reg,comm) - p_marketData(reg,"qDom",comm,"cal");

*   - For supply quantities
p_problemDiff(reg,"na","qSup",comm,"cal")
    = v_qSup.L(reg,comm) - p_marketData(reg,"qSup",comm,"cal");

*   - For demand prices
p_problemDiff(reg,"na","pConsDom",comm,"cal")
    = v_pConsDom.L(reg,comm) - p_marketData(reg,"pConsDom",comm,"cal");

*   - For traded quantities
p_problemDiff(reg,reg1,"qTrd",comm,"cal")
    = v_qTrd.L(reg,reg1,comm) - p_tradeData(reg,reg1,"qTrd",comm,"cal");

*   - Utility of imports match calibration point
p_problemDiff(reg,"na","uImp",comm,"cal")
    = p_marketData(reg,"uImp",comm,"cal") - v_uImp.L(reg,comm);

*   - Price index of imports match calibrated value
p_problemDiff(reg,"na","pImpAgg",comm,"cal")
    = v_armImpPriceIndex.L(reg,comm) - p_marketData(reg,"pImpAgg",comm,"cal");

problemData(reg,regtot,cols,comm,"cal")
    $ [abs(p_problemDiff(reg,regtot,cols,comm,"cal")) gt p_problemTol] = yes;

$batinclude "assert_that_set_is_empty.gms" problemData "Some price or quantity is not calibrating correctly" %ERROR_FILE%


*-------------------------------------------------------------------------------
* TESTS for plausible simulation behaviour
* "Only code if there is a failing test"
*-------------------------------------------------------------------------------

* --- Increase total demand (exogenous) in all regions vs calibration and verify
*     plausibility of results.
option kill=problemData;
option kill=p_problemDiff;
option kill=colsToCheck;

colsToCheck("qSup") = yes;
colsToCheck("pProd") = yes;
colsToCheck("pConsDom") = yes;


p_marketData(reg,"uDem",comm,"sim") = p_marketData(reg,"uDem",comm,"cal")*1.1;
p_marketData(reg,"ppMrg",comm,"sim") = p_marketData(reg,"ppMrg",comm,"cal")*1.0;
solve m_trade using cns;

*   - Total supply should increase in all regions vs calibration
p_problemDiff(reg,"na","qSup",comm,"sim")
    = v_qSup.L(reg,comm) - p_marketData(reg,"qSup",comm,"cal");

*   - Producer price should increase in all regions vs calibration
p_problemDiff(reg,"na","pProd",comm,"sim")
    = v_pProd.L(reg,comm) - p_marketData(reg,"pProd",comm,"cal");

*   - The consumer price of domestic goods should increase in all regions vs calibration
p_problemDiff(reg,"na","pConsDom",comm,"sim")
    = v_pConsDom.L(reg,comm) - p_marketData(reg,"pConsDom",comm,"cal");

*   - The export price should increase in all regions vs calibration
p_problemDiff(reg,"na","pFob",comm,"sim")
    = v_pFob.L(reg,comm) - p_marketData(reg,"pFob",comm,"cal");

*   - The import price (pCif) should increase in all regions vs calibration
p_problemDiff(reg,"na","pFob",comm,"sim")
    = v_pFob.L(reg,comm) - p_marketData(reg,"pFob",comm,"cal");


problemData(reg,"na",colsToCheck,comm,"sim")
    $ [p_problemDiff(reg,"na",colsToCheck,comm,"sim") lt p_problemTol] = yes;



$batinclude "assert_that_set_is_empty.gms" problemData "Behavioural problem: increased demand does not provoke increased supply." %ERROR_FILE%
