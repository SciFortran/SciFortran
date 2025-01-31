#!/bin/bash

# Function to convert RGB to HEX
rgb_to_hex() {
    printf "#%02X%02X%02X" "$1" "$2" "$3"
}

# Read the file line by line
while IFS= read -r line; do
    # Check if the line matches the pattern "type(rgb_color),parameter :: green = rgb_color(rr, gg, bb)"
    if [[ "$line" =~ rgb_color\([[:space:]]*([0-9]+)[[:space:]]*,[[:space:]]*([0-9]+)[[:space:]]*,[[:space:]]*([0-9]+)\) ]]; then
        # Extract the RGB values (rr, gg, bb)
        rr="${BASH_REMATCH[1]}"
        gg="${BASH_REMATCH[2]}"
        bb="${BASH_REMATCH[3]}"
        
        # Convert the RGB values to HEX
        hex_color=$(rgb_to_hex "$rr" "$gg" "$bb")
        
        # Append the HTML comment with the HEX color
        echo "$line !:raw-html:\`<span class=\"colored-square\" style=\"background-color: $hex_color;\"></span>\`"
    else
        # If the line doesn't match the pattern, just print it as is
        echo "$line"
    fi
done < "SF_COLORS.f90" > "SF_COLORS_NEW.f90"
