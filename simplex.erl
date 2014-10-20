-module(simplex).

-export([noise/3, setup/0]).

-define(F3, 1.0/3.0).
-define(G3, 1.0/6.0).


setup() ->
        P = [151,160,137,91,90,15,
  131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
  190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
  88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
  77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
  102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
  135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
  5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
  223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
  129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
  251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
  49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
  138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180],

        Perm = [ lists:nth((Per band 255) + 1, P) || Per <- lists:seq(0, 511)],
        PermMod12 = [Per12 rem 12 || Per12 <- Perm],
	put(perm,Perm),
	put(permMod, PermMod12),
	ok.

dot({GX, GY, GZ}, X, Y, Z) ->
	GX*X + GY*Y + GZ*Z.

floor(X) ->
    T = erlang:trunc(X),
    case (X - T) of
        Neg when Neg < 0 -> T - 1;
        Pos when Pos > 0 -> T;
        _ -> T
    end.

noise(Xin, Yin, Zin) ->
	% Setup

	Grad = [
		{1,1,0}, {-1,1,0}, {1,-1,0}, {-1,-1,0},
		{1,0,1}, {-1,0,1}, {1,0,-1}, {-1,0,-1},
		{0,1,1}, {0,-1,1}, {0,1,-1}, {0,-1,-1}
	],
	Perm = get(perm),
	PermMod12 = get(permMod),
	% Actual Code	
	S = (Xin+Yin+Zin)*?F3,

	I = floor(Xin + S),
	J = floor(Yin + S),
	K = floor(Zin + S),

	T = (I+J+K) * ?G3,

	X = I - T,
	Y = J - T,
	Z = K - T,

	X0 = Xin - X,
	Y0 = Yin - Y,
	Z0 = Zin - Z,

	{I1, J1, K1, I2, J2, K2} = if
		X0 >= Y0 ->
			if
				Y0 >= Z0 -> {1,0,0, 1,1,0};
				X0 >= Z0 -> {1,0,0, 1,0,1};
				true -> {0,0,1, 1,0,1}
			end;
		true ->
			if
				Y0 < Z0 -> {1,0,0, 1,1,0};
				X0 < Z0 -> {0,1,0, 0,1,1};
				true -> {0,1,0, 1,1,0}
			end
	end,

	X1 = X0 - I1 + ?G3,
	Y1 = Y0 - J1 + ?G3,
	Z1 = Z0 - K1 + ?G3,
	
	X2 = X0 - I1 + 2.0 * ?G3,
	Y2 = Y0 - J1 + 2.0 * ?G3,
	Z2 = Z0 - K1 + 2.0 * ?G3,

	X3 = X0 - 1.0 + 3.0 * ?G3,
	Y3 = Y0 - 1.0 + 3.0 * ?G3,
	Z3 = Z0 - 1.0 + 3.0 * ?G3,

	II = erlang:trunc(I) band 255,
	JJ = erlang:trunc(J) band 255,
	KK = erlang:trunc(K) band 255,

	P0 = getPerm(KK, Perm),
	P1 = getPerm(JJ + P0, Perm),
	GI0 = getPerm(II + P1, PermMod12),

	P2 = getPerm(KK + K1, Perm),
	P3 = getPerm(JJ + J1 + P2, Perm),
	GI1 = getPerm(II + I1 + P3, PermMod12),

	P4 = getPerm(KK + K2, Perm),

	P5 = getPerm(JJ + J2 + P4, Perm),
	GI2 = getPerm(II + I2 + P5, PermMod12),

	P6 = getPerm(KK + 1, Perm),
	P7 = getPerm(JJ + 1 + P6, Perm),
	GI3 = getPerm(II + 1 + P7, PermMod12),

	T0 = 0.6 - X0*X0 - Y0*Y0 - Z0*Z0,

	N0 = if
		T0 < 0 -> 0.0;
		true -> 
			T00 = T0*T0,
			T00*T00 * dot(lists:nth(GI0 + 1, Grad), X0, Y0, Z0)
	end,

	T1 = 0.6 - X1*X1 - Y1*Y1 - Z1*Z1,

	N1 = if
		T1 < 0 -> 0.0;
		true ->
			T11 = T1*T1,
			T11*T11 * dot(lists:nth(GI1 + 1, Grad), X1, Y1, Z1)
	end,

	T2 = 0.6 - X2*X2 - Y2*Y2 - Z2*Z2,

        N2 = if
                T2 < 0 -> 0.0;
                true -> 
                        T22 = T2*T2,
                        T22*T22 * dot(lists:nth(GI2 + 1, Grad), X2, Y2, Z2)
        end,

	T3 = 0.6 - X3*X3 - Y3*Y3 - Z3*Z3,

        N3 = if
                T3 < 0 -> 0.0;
                true -> 
                        T33 = T3*T3,
                        T33*T33 * dot(lists:nth(GI3 + 1, Grad), X3, Y3, Z3)
        end,

	32.0 * (N0+N1+N2+N3).

getPerm(0, [H|_]) -> H;
getPerm(Left,[]) -> 
	Overflow = get(max_overflow),
	if
		Overflow == undefined orelse Overflow < Left -> put(max_overflow,Left);
		true -> ok
	end,
	0;
getPerm(Index, [_|T]) -> 
	getPerm(Index - 1, T).
