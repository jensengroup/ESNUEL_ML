# A script for reference how the carb_acids dataset
# was generated based on the paper: https://onlinelibrary.wiley.com/doi/full/10.1002/ffj.3327
# Haven't tested it (imports are missing etc.).
# Needs my small molecule package dependencies to work, too.
# But just so we know what we did to the data.
#
# It is: in days, and: with np.log10 applied to that.


if __name__ == "__main__":
  csv = StringIO(
  """
  Compound	tAA	tA	tN	tB
  (3H)-Ethylfuranone	110	160	93	2.2
  (5H)-Ethylfuranone	44,000 (120 y)	660,000 (1,800 y)	40,000 (110 y)	400
  (E,E)-Farnesyl acetate	18	290	91	0.91 (22 h)
  (E)-Oak lactone	1,100	16,000 (44 y)	9,500	97
  (E)-Whiskey lactone	940	13,000 (36 y)	7,300	79
  (Z)-3-Hexenyl hexanoate	45	700	820	8.3
  (Z)-6-Dodecene-γ-lactone	560	8,700	8,100	83
  (Z)-Whiskey lactone	940	13,000 (36 y)	7,300	79
  2-Methylbutyl acetate	130	2,000	2,500	26
  3-Mercaptohexyl acetate	99	1,500	1,300	13
  4-Carbethoxybutyrolactone (chaina)	260	2,300	200	2.1
  4-Carbethoxybutyrolactone (ringb)	550	8,400	2,700	27
  4-Hexanolide	980	15,000 (41 y)	8,800	90
  5-Octanolide	1,600	24,000 (66 y)	14,000 (38 y)	140
  7-Methoxycoumarin	43,000 (120 y)	310,000 (850 y)	1,100	11
  Benzyl acetate	42	660	160	1.6
  Benzyl benzoate	17,000 (47 y)	200,000 (550 y)	770	7.7
  Benzyl butanoate	130	2,100	400	4.1
  Bornyl benzoate	130,000 (360 y)	1,800,000 (4,900 y)	17,000 (47 y)	170
  Bornyl butyrate	1,000	16,000 (44 y)	9,200	94
  Bornyl formate	22	56	14	0.18 (4.3 h)
  Bornyl isovalerate	1,500	23,000 (63 y)	13,000 (36 y)	140
  Butyl acetate	49	780	1,100	11
  Butyl benzoate	21,000 (58 y)	300,000 (820 y)	5,400	54
  Butyl decanoate	190	3,000	3,300	34
  Butyl hexanoate	190	2,900	3,200	33
  Butyl laurate	190	3,000	3,300	34
  Butyl methylbutyrate	720	9,500	8,600	100
  Butyl octanoate	190	2,900	3,200	34
  Carvyl acetate	67	1,100	300	3.1
  Citronellyl acetate	100	1,600	2,000	20
  Citronellyl butyrate	320	4,900	5,100	53
  Citronellyl isobutyrate	500	6,900	6,400	74
  Citronellyl valerate	360	5,500	5,700	60
  Coumarin	28,000 (77 y)	150,000 (410 y)	300	3.2
  Diethyl 2-hydroxyglutarate	860	3,600	520	5.8
  Diethyl malate	180	2,700	580	5.9
  Diethyl malonate	30	470	26	0.27 (6.5 h)
  Diethyl succinate	36	570	260	2.6
  Diethyl tartrate	1,200	2,200	320	3.7
  Dihydrocarvyl acetate	200	3,200	2,200	22
  Ethyl (E)-cinnamate	2,400	36,000 (99 y)	1,000	10
  Ethyl 3-hydroxybutanoate	190	2,800	1,700	18
  Ethyl 3-methylbutanoate	200	3,100	3,400	35
  Ethyl acetate	42	660	930	9.4
  Ethyl benzoate	18,000 (49 y)	260,000 (700 y)	4,800	48
  Ethyl butanoate	140	2,100	2,400	25
  Ethyl cinnamate	2,400	36,000 (99 y)	1,000	10
  Ethyl cyclohexanoate	770	9,300	4,300	48
  Ethyl decanoate	160	2,500	2,800	29
  Ethyl dihydrocinnamate	190	2,900	2,000	21
  Ethyl formate	11	26	11	0.17 (4.1 h)
  Ethyl heptylate	160	2,500	2,800	29
  Ethyl hexadecanoate	160	2,500	2,800	29
  Ethyl hexanoate	160	2,500	2,800	29
  Ethyl hydroxybutanoate	160	2,400	2,200	22
  Ethyl hydroxyhexanoate	760	11,000 (30 y)	5,900	62
  Ethyl isobutyrate	220	3,100	3,100	34
  Ethyl isohexanoate	140	2,200	2,500	25
  Ethyl lactate	240	1,600	230	2.5
  Ethyl laurate	160	2,500	2,800	29
  Ethyl mercaptopropionate	280	1,400	170	1.9
  Ethyl methylbutyrate	620	8,300	7,500	88
  Ethyl octanoate	160	2,500	2,800	29
  Ethyl phenylacetate	180	2,800	610	6.1
  Ethyl propionate	96	1,500	1,800	18
  Ethyl salicylate	65,000 (180 y)	870,000 (2,400 y)	10,000 (27 y)	100
  Ethyl tetradecanoate	160	2,500	2,800	29
  Ethyl undecanoate	160	2,500	2,800	29
  Ethyl valerate	150	2,400	2,700	27
  Ethyl vanillate	41,000 (110 y)	620,000 (1,700 y)	29,000 (79 y)	290
  Geranyl acetate	18	290	91	0.91 (22 h)
  Geranyl butyrate	43	670	240	2.4
  Geranyl isovalerate	64	1,000	340	3.4
  Geranyl valerate	48	760	270	2.7
  Hexenyl acetate	30	470	270	2.7
  Hexyl acetate	50	800	1,100	11
  Hexyl butanoate	160	2,500	2,800	29
  Hexyl hexanoate	190	2,900	3,200	34
  Hexyl methylbutyrate	740	9,700	8,700	100
  Hexyl octanoate	200	3,000	3,300	34
  Isoamyl acetate	35	550	780	7.9
  Isobornyl formate	22	56	14	0.18 (4.3 h)
  Isobornyl propionate	730	11,000 (30 y)	6,700	69
  Isobutyl acetate	53	830	1,100	11
  Isopropyl benzoate	34,000 (93 y)	510,000 (1,400 y)	8,400	84
  Isopropyl hexanoate	310	4,900	5,100	52
  Isopropyl palmitate	320	5,000	5,200	53
  Isopulegyl acetate	210	3,300	2,200	22
  Linalyl acetate	420	6,600	1,300	13
  Linalyl butyrate	1,400	21,000 (58 y)	3,500	35
  Linalyl formate	11	120	1.8	0.018 (26 min)
  Linalyl isovalerate	2,000	32,000 (88 y)	5,000	50
  Linalyl valerate	1,500	24,000 (66 y)	3,900	39
  Mercaptomethylbutyl formate	9.7	18	7.2	0.11 (2.6 h)
  Methyl 2-(methylthio)acetate	130	660	38	0.40 (9.6 h)
  Methyl 2-methylbutanoate	410	5,000	4,800	62
  Methyl 2-methylpropanoate	140	2,000	2,000	24
  Methyl 3-methylbutanoate	140	2,000	2,400	25
  Methyl anthranilate	63,000 (170 y)	970,000 (2,700 y)	37,000 (100 y)	370
  Methyl benzoate	12,000 (33 y)	180,000 (490 y)	3,500	35
  Methyl butanoate	90	1,400	1,600	17
  Methyl cinnamate	1,600	24,000 (66 y)	700	7.0
  Methyl cyclohexanecarboxylate	510	5,400	2,800	34
  Methyl decanoate	110	1,700	1,900	20
  Methyl dihydroepijasmonate	160	2,400	1,900	20
  Methyl epijasmonate	150	2,300	1,700	18
  Methyl geranate	9,200	130,000 (360 y)	71,000 (190 y)	750
  Methyl hexanoate	110	1,600	1,900	20
  Methyl jasmonate	150	2,300	1,700	18
  Methyl laurate	110	1,700	1,900	20
  Methyl nonanoate	110	1,700	1,900	20
  Methyl octadecenoate	110	1,700	1,900	20
  Methyl octanoate	110	1,700	1,900	20
  Methyl salicylate	44,000 (120 y)	590,000 (1,600 y)	7,700	77
  Methyl tetradecanoate	110	1,700	1,900	20
  Methyl vanillate	28,000 (77 y)	410,000 (1,100 y)	21,000 (58 y)	210
  Neryl acetate	18	290	91	0.91 (22 h)
  Nonyl acetate	50	800	1,100	11
  Octyl acetate	420	6,500	7,300	74
  p-Menth-1-en-9-yl acetate	300	4,700	5,500	56
  Pantolactone	2,300	6,300	1,000	12
  Pentyl butanoate	160	2,500	2,800	29
  Phenylethyl benzoate	20,000 (55 y)	290,000 (790 y)	4,100	41
  Propyl butyrate	150	2,400	2,700	27
  Propyl hexanoate	180	2,800	3,100	32
  Propyl propanoate	110	1,700	1,900	20
  R-δ-Decenolactone	80,000 (220 y)	1,200,000 (3,300 y)	89,000 (240 y)	890
  Sotolon	130,000 (360 y)	1,200,000 (3,300 y)	170,000 (470 y)	1,900
  Terpinyl acetate	3,700	58,000 (160 y)	33,000 (90 y)	330
  Wine lactone	5,100	58,000 (160 y)	22,000 (60 y)	250
  β-Phenethyl acetate	47	750	820	8.3
  γ-Butyrolactone	540	8,100	8,100	85
  γ-Decalactone	980	15,000 (41 y)	8,800	90
  γ-Dodecalactone	980	15,000 (41 y)	8,800	90
  γ-Nonalactone	980	15,000 (41 y)	8,800	90
  γ-Octalactone	980	15,000 (41 y)	8,800	90
  γ-Undecalactone	980	15,000 (41 y)	8,800	90
  δ-Decalactone	1,600	24,000 (66 y)	14,000 (38 y)	140
  δ-Dodecalactone	1,600	24,000 (66 y)	14,000 (38 y)	140
  δ-Undecalactone	1,600	24,000 (66 y)	14,000 (38 y)	140
  """
  )
  df = pd.read_csv(csv,sep="\t")
  
  def time_parse(val):
      return float(val.split()[0].replace(",",""))
  
  for col in ["tAA","tA","tN","tB"]:
      df[col] = df[col].apply(time_parse)
  
  df["canonical_smiles"] = df["Compound"].apply(iupac_to_smiles)
  
  print(len(df))
  df["canonical_smiles"] = df["canonical_smiles"].apply(lambda s: s.strip())
  df = df[df.canonical_smiles.apply(len) > 0]
  print(len(df))
  
  df.to_csv("/home/gnlop/carb_acids/carb_acids.csv")
  df = pd.read_csv(data / "carb_acids.csv",sep="\t")
  df["canonical_smiles"] = df["Compound"].apply(iupac_to_smiles)
  df["canonical_smiles"] = df["canonical_smiles"].apply(lambda s: s.strip())
  df = df[df.canonical_smiles.apply(len) > 0]
  add_cv_by_col(df,"canonical_smiles")
  df["split"] = df["cv"]
  def time_parse(val):
      return float(val.split()[0].replace(",",""))

  for col in ["tAA","tA","tN","tB"]:
      df[col] = df[col].apply(time_parse)
      df[col] = df[col].apply(np.log10)
  df.to_csv(data / "carb_acids_2.csv")
