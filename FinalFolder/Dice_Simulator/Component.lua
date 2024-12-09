
-- Function to adjust the acceleration
function AdjustAcceleration(x, y, z)
    updateAcceleration(x, y, z)  -- Calls the C++ function
    print(string.format("Acceleration updated to: (%f, %f, %f)", x, y, z))
end

-- Edit these values to edit the particles acceleration
AdjustAcceleration(0.0, -1.0, 1.0)
