let n = 243
let a = 6
let x = 7
let z = 15

let q = ((x + z) * (x + z)) mod n
let r = ((x - z) * (x - z)) mod n
let s = (q - r) mod n

let _ = print_int q;;
let _ = print_string "\n";;
let _ = print_int r;;
let _ = print_string "\n";;
let _ = print_int s;;
let _ = print_string "\n";;

let xm = (q * r) mod n
let zm = (s * (r + (a + 2) * s / 4)) mod n

let _ = print_int xm;;
let _ = print_string "\n";;
let _ = print_int zm;;
let _ = print_string "\n";;
