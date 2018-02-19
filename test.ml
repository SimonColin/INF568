let n = 101
let a = 38
let x = 2
let z = 1

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
