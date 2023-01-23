#pragma once

std::random_device rd;
std::mt19937 rng(rd());
std::uniform_real_distribution<> dist(0, 1);
std::uniform_int_distribution<> distInt(1, n + 1); // define the range
std::uniform_int_distribution<> distInt01(0, 1); // define the range