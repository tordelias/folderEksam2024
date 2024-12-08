-- Component properties
ComponentProperties = {
    acceleration = { x = 0.0, y = -9.8, z = 0.0 }
}

-- Function to adjust the acceleration and lifetime
function AdjustAcceleration(x, y, z)
    updateAcceleration(x, y, z)  -- Calls the C++ function
    print(string.format("Acceleration updated to: (%f, %f, %f) with lifetime: %f", x, y, z))
end

-- Edit these values to edit the particles acceleration and lifetime
AdjustAcceleration(0.0, -1.0, 1.0)
