module ADF_test =
    
    open System
    open MathNet.Numerics.LinearRegression
    open MathNet.Numerics.LinearAlgebra

    type RegressionType =
        | Constant
        | ConstantTrend
        | ConstantLinearQuadraticTrend
        | NoConstantNoTrend


    module private PValue =

        open MathNet.Numerics.Distributions

        let private  tau_star_nc = [| -1.04; -1.53; -2.68; -3.09; -3.07; -3.77 |]
        let private  tau_min_nc = [| -19.04;-19.62;-21.21;-23.25;-21.63;-25.74 |]
        let private  tau_max_nc = [| Double.PositiveInfinity;1.51;0.86;0.88;1.05;1.24 |]
        let private  tau_star_c = [| -1.61; -2.62; -3.13; -3.47; -3.78; -3.93 |]
        let private  tau_min_c = [| -18.83;-18.86;-23.48;-28.07;-25.96;-23.27 |]
        let private  tau_max_c = [| 2.74;0.92;0.55;0.61;0.79;1.0 |]
        let private  tau_star_ct = [| -2.89; -3.19; -3.50; -3.65; -3.80; -4.36 |]
        let private  tau_min_ct = [| -16.18;-21.15;-25.37;-26.63;-26.53;-26.18 |]
        let private  tau_max_ct = [| 0.7;0.63;0.71;0.93;1.19;1.42 |]
        let private  tau_star_ctt = [| -3.21;-3.51;-3.81;-3.83;-4.12;-4.63 |]
        let private  tau_min_ctt = [| -17.17;-21.1;-24.33;-24.03;-24.33;-28.22 |]
        let private  tau_max_ctt = [| 0.54;0.79;1.08;1.43;3.49;1.92 |]
        let private  tau_nc_smallp =  [| 0.6344; 1.2378; 0.032496 |]
        let private  tau_c_smallp = [| 2.1659; 1.4412; 0.038269 |]
        let private  tau_ct_smallp = [| 3.2512; 1.6047; 0.049588 |]
        let private  tau_ctt_smallp = [| 4.0003; 1.658; 0.048288 |]
        let private  tau_nc_largep = [| 0.4797; 0.93557; -0.06999; 0.033066 |]
        let private  tau_c_largep = [| 1.7339; 0.93202; -0.12745; -0.010368 |]
        let private  tau_ct_largep = [| 2.5261; 0.61654; -0.37956; -0.060285 |]
        let private  tau_ctt_largep = [| 3.0778; 0.49529; -0.41477; -0.059359 |]


        let getTStatPValue regType tstat =
    
            let maxStat, minStat, starStat, smallp, largep =
                match regType with
                | NoConstantNoTrend -> 
                    (tau_max_nc, tau_min_nc, tau_star_nc, tau_nc_smallp, tau_nc_largep)
                | Constant -> 
                    (tau_max_c, tau_min_c, tau_star_c, tau_c_smallp, tau_c_largep)
                | ConstantTrend -> 
                    (tau_max_ct, tau_min_ct, tau_star_ct, tau_ct_smallp, tau_ct_largep)
                | ConstantLinearQuadraticTrend -> 
                    (tau_max_ctt, tau_min_ctt, tau_star_ctt, tau_ctt_smallp, tau_ctt_largep)

            if tstat > maxStat.[0] then 1.0
            elif tstat < minStat.[0] then 0.0
            else
        
                let tau_coeff =
                    match tstat <= starStat.[0] with
                    | true -> smallp
                    | false -> largep

                let polyval =
                    tau_coeff
                    |> Array.mapi(fun i c -> c * (tstat ** (i |> float)))
                    |> Array.sum

                Normal.CDF(-polyval, 1.0, 0.0)




    type MaxLags = MaxLags of int


    type LagCalculation =
        | Default
        | AIC
        | BIC
        | TStat


    type ADFparameters =
        {
            RegressionType: RegressionType
            LagCalculation: LagCalculation
            MaxLags: MaxLags option

        }


    type private RegressionResults =
        {
            Coefficient: float
            StdError: float
            TStat: float
        }


    type private SSEresult =
        {
            SSE: float
            k: int
        }


    // Given an endogenous vector and matrix of exogenous variables
    // run a multiple regression to fit a linear model
    // and return the regression stats
    let private getRegressionResults hasIntercept (y:Vector<float>) (x:Matrix<float>)  =
       
        let coeffs = 
            MultipleRegression.DirectMethod(x.ToRowArrays(), y.ToArray(), hasIntercept, DirectRegressionMethod.QR)
            |> DenseVector.ofArray

        let resids = (y - (x * coeffs))

        let sigmaSquared = (resids .* resids |> Vector.sum) / ((x.RowCount - x.ColumnCount) |> float)

        let inv = Matrix.transpose(x) * x |> Matrix.inverse

        let mCovar = sigmaSquared * inv

        let stdErrors = mCovar.Diagonal().ToArray() |> Array.map(fun v -> v ** 0.5)

        (coeffs.ToArray(), stdErrors)
        ||> Array.zip
        |> Array.map(fun (coeff, stdError) ->
            
            let tstat = coeff / stdError

            {
                Coefficient = coeff
                StdError = stdError
                TStat = tstat
            }
            
            )


   
    // for a given max lag, check that it is in bounds vs the number of observations
    // and for no max lag, calculate the max lag
    let private calcMaxLag lagCalc regType nObs maxLags =
        
        let ntrend =
            match regType with
            | Constant -> 1
            | ConstantTrend -> 2
            | ConstantLinearQuadraticTrend -> 3
            | NoConstantNoTrend -> 0

        match lagCalc, maxLags with
        | _, Some (MaxLags lag) ->
            if lag > ((nObs / 2) - ntrend - 1) then
                failwith ""
            else
                lag
        | _ ->
            (12.0 * ((nObs |> float) / 100.0) ** 0.25)
            |> Math.Ceiling
            |> int
            |> min ((nObs / 2) - ntrend - 1)
            


    // Given an endogenous vector and matrix of exogenous variables
    // run a multiple regression to fit a linear model
    // and return the sum of the squared residuals and the number of coefficients
    let private calcSSE (y:Vector<float>) (x:Matrix<float>) =
        
        let coeffs = 
            MultipleRegression.DirectMethod(x, y, DirectRegressionMethod.QR)

        let resids = (y - (x * coeffs))

        let sse = resids.ToArray() |> Array.map(fun r -> r * r) |> Array.sum

        { SSE = sse; k = coeffs.Count }


    // Given an endogenous vector and matrix of exogenous variables
    // calculate the best number of lags to use based on the lag calculation
    let private calcBestLag (x:Matrix<float>) (y:Vector<float>) startLag maxLag lagCalc =
        
        match lagCalc with
        | Default -> maxLag
        | TStat ->

            let fit =
                [| for l in (startLag + maxLag)..(-1)..startLag do
                    
                    let regressionResults =
                        x.GetSlice(None,None,None, l |> Some)
                        |> getRegressionResults false y

                    yield (l - startLag , regressionResults.[l].TStat)
                |]

            match fit |> Array.tryFind(fun (i,fit) -> Math.Abs(fit) > 1.6448536269514722) with
            | Some (i,fit) -> i
            | None -> fit |> Array.maxBy(fun (i,fit) -> Math.Abs(fit)) |> fst

        | AIC | BIC ->
            
            let fit =
                [| for l in startLag..(startLag + maxLag) do
                    
                    yield
                        x.GetSlice(None, None, None, l |> Some)
                        |> calcSSE y
                    |]

            let n = y.Count |> float

            match lagCalc with
            | AIC ->

                fit
                |> Array.mapi(fun i f -> (i,f))
                |> Array.minBy(fun (nlags, f) ->
                        
                    let k = f.k |> float
                    let sse = f.SSE

                    let aic = n * (Math.Log(sse / n)) + (2.0 * k)
                    aic)
                |> fst

            | BIC ->
            
                fit
                |> Array.mapi(fun i f -> (i,f))
                |> Array.minBy(fun (i,f) ->
                        
                    let k = f.k |> float
                    let sse = f.SSE

                    let bic = n * (Math.Log(sse / n)) + k * Math.Log(n) 
                    bic)
                |> fst

            | _ -> failwith "Should not get here"


    let private addTrend regression x =
        
            match regression with
            | NoConstantNoTrend -> x
            | Constant ->
                x
                |> Array.map(fun terms -> 
                    terms |> Array.append [| 1.0 |])
            | ConstantTrend ->
                x
                |> Array.mapi(fun i terms -> 
                    terms |> Array.append [| 1.0; (i + 1 |> float) |])
            | ConstantLinearQuadraticTrend ->
                x
                |> Array.mapi(fun i terms -> 
                    terms 
                    |> Array.append [| 1.0; (i + 1 |> float); (i + 1 |> float) ** 2.0 |])


    let private lagMatrix (x:float[]) nLags =
    

        let xDiff =
            x |> Array.pairwise |> Array.map(fun (x1, x2) -> x2 - x1)

        match nLags with
        | 0 -> 
            let xdAll = 
                x 
                |> Array.take (x.Length - 1)
                |> Array.map(fun x -> [| x |])
        
            (xdAll, xDiff)
        | _ ->
        
            let lagTerms =
                    xDiff.[nLags - 1..]
                    |> Array.mapi(fun i _ -> xDiff.[i..(i + nLags - 1)] |> Array.rev )

            let xdAll =
                (x.[nLags..], lagTerms) 
                ||> Array.map2(fun x lags -> lags |> Array.append [| x |] )
                |> Array.take (lagTerms.Length - 1)

            let xdShort = xDiff.[nLags..]

            (xdAll, xdShort)


    let private calcTstat (xdAll:float[][]) xdShort regression maxLag =
    
                let y = xdShort |> DenseVector.ofArray
                let x = 
                    xdAll 
                    |> Array.map(fun terms -> terms.[0..maxLag])
                    |> addTrend regression
                    |> Array.transpose
                    |> DenseMatrix.ofColumnArrays

                let regressionResults =
                    getRegressionResults false y x

                let i =
                    match regression with
                    | NoConstantNoTrend -> 0
                    | Constant -> 1
                    | ConstantTrend -> 2
                    | ConstantLinearQuadraticTrend -> 3

                regressionResults.[i].TStat


    let ADFTest (adfParams:ADFparameters) (data:float[]) =
        
        let lagCalc = adfParams.LagCalculation
        let regType = adfParams.RegressionType
        let maxLags = adfParams.MaxLags

        let nObs = data.Length

        let maxLag = calcMaxLag lagCalc regType nObs maxLags

        // create initial lag matrix and endogenous vector
        let xdAll, xdShort = lagMatrix data maxLag

        let tstat =
            
            match lagCalc with
            | Default ->

                // if we are simply using the maxLag to determine the model
                // then use the initial lag matrix and endogenous vector
                calcTstat xdAll xdShort regType maxLag

            | lagCalc ->
            
                // add constant and trend components as required
                let fullRHS = addTrend regType xdAll

                let startLag =
                    match regType with
                    | Constant -> 1
                    | ConstantTrend -> 2
                    | ConstantLinearQuadraticTrend -> 3
                    | NoConstantNoTrend -> 0

                // determine the best lag based on the lag calc
                let bestLag = 
                    calcBestLag 
                        (fullRHS |> DenseMatrix.ofRowArrays)
                        (xdShort |> DenseVector.ofArray)
                        startLag maxLag lagCalc

                // build updated lag matrix and endogenous vector
                // using the best lag
                let xdAll, xdShort = lagMatrix data bestLag
            
                calcTstat xdAll xdShort regType bestLag

        PValue.getTStatPValue regType tstat
