(* independent of lsvir.c but does the same thing *)

ring3 = {3,
  4 - 3*Sqrt[3]/Pi,
  15/8,
  4 - 9*Sqrt[3]/2/Pi,
  159/2^7,
  4 - 27*Sqrt[3]/5/Pi,
  3*289/2^10,
  4 - 3*279*Sqrt[3]/140/Pi,
  3*6413/2^15,
  4 - 3*297*Sqrt[3]/140/Pi,
  3*35995/2^18,
  4 - 3*243*Sqrt[3]/110/Pi
};

ring = {16/3,
  8 - 128/3/Pi^2,
  272/105,
  8 - 8192/135/Pi^2,
  4016/3003,
  8 - 65536/945/Pi^2,
  296272/415701,
  8 - 134217728/1819125/Pi^2,
  1234448/3187041,
  8 - 1744830464/22920975/Pi^2,
  55565456/260468169,
  8 - 4312147165184/55720890225/Pi^2
};

diam = {-14/3,
  -8 + 8*Sqrt[3]/Pi + 20/3/Pi^2,
  -6347/3360,
  -8 + 12*Sqrt[3]/Pi + 173/135/Pi^2,
  -20830913/24600576,
  -8 + 72*Sqrt[3]/5/Pi - 193229/37800/Pi^2,
  -87059799799/217947045888,
  -8 + 558*Sqrt[3]/35/Pi - 76667881/7276500/Pi^2,
  -332647803264707/1711029608251392,
  -8 + 594*Sqrt[3]/35/Pi - 9653909/654885/Pi^2,
  -865035021570458459/8949618140032008192,
  -8 + 972*Sqrt[3]/55/Pi - 182221984415/10188962784/Pi^2
};

star = {4,
  8 - 12*Sqrt[3]/Pi + 8/Pi^2,
  (-89/280 - 219*Sqrt[2]/1120/Pi + 4131/2240ArcCos[1/3]/Pi)*4,
  8 - 18*Sqrt[3]/Pi + 238/9/Pi^2,
  (-163547/128128 - 3888425*Sqrt[2]/8200192/Pi + 67183425/16400384ArcCos[1/3]/Pi)*4,
  8 - 108*Sqrt[3]/5/Pi + 37259/900/Pi^2,
  (-283003297/141892608 - 159966456685*Sqrt[2]/217947045888/Pi + 292926667005/48432676864ArcCos[1/3]/Pi)*4,
  8 - 837*Sqrt[3]/35/Pi + 5765723/110250/Pi^2,
  (-88041062201/34810986496 - 2698457589952103*Sqrt[2]/2851716013752320/Pi + 8656066770083523/1140686405500928ArcCos[1/3]/Pi)*4,
  8 - 891*Sqrt[3]/35/Pi + 41696314/694575/Pi^2,
  (-66555106087399/22760055898112 - 16554115383300832799*Sqrt[2]/14916030233386680320/Pi + 52251492946866520923/5966412093354672128ArcCos[1/3]/Pi)*4,
  8 - 1458*Sqrt[3]/55/Pi + 88060381669/1344697200/Pi^2
};

For [ dim = 1, dim <= 12, dim++,
  Z = Apart[ (ring[[dim]]*3 - star[[dim]]*2) / 8 ];
  Print["D ", dim, ", Z4/V^3 ", (Z // InputForm), " == ", N[Z, 20]];
];

For [ dim = 1, dim <= 12, dim++,
  B = Apart[ -3*(ring[[dim]]/8 + diam[[dim]]/4 + star[[dim]]/24) ];
  Print["D ", dim, ", B4/B2^3 ", (B // InputForm), " == ", N[B, 20]];
];
